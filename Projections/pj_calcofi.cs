// Conversions for the California Cooperative Oceanic Fisheries Investigations
// Line/Station coordinate system following the algorithm of:
// Eber, L.E., and  R.P. Hewitt. 1979. Conversion algorithms for the CALCOFI
// station grid. California Cooperative Oceanic Fisheries Investigations Reports
// 20:135-137. (corrected for typographical errors).
// http://www.calcofi.org/publications/calcofireports/v20/Vol_20_Eber___Hewitt.pdf
//
// They assume 1 unit of CalCOFI Line == 1/5 degree in longitude or
// meridional units at reference point O, and similarly 1 unit of CalCOFI
// Station == 1/15 of a degree at O.
//
// By convention, CalCOFI Line/Station conversions use Clarke 1866 but we use
// whatever ellipsoid is provided.

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_calcofi : PJ
	{
		public override string Name { get { return "calcofi"; } }
		public override string DescriptionName { get { return "Cal Coop Ocean Fish Invest Lines/Stations"; } }
		public override string DescriptionType { get { return "Cyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double EPS10=1.0e-10;
		const double DEG_TO_LINE=5;
		const double DEG_TO_STATION=15;
		const double LINE_TO_RAD=0.0034906585039886592;
		const double STATION_TO_RAD=0.0011635528346628863;
		const double PT_O_LINE=80; // reference point O is at line 80,
		const double PT_O_STATION=60; // station 60,
		const double PT_O_LAMBDA=-2.1144663887911301; // lon -121.15 and
		const double PT_O_PHI=0.59602993955606354; // lat 34.15
		const double ROTATION_ANGLE=0.52359877559829882; // CalCOFI angle of 30 deg in rad

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double oy; // pt O y value in Mercator
			double l1, l2; // l1 and l2 are distances calculated using trig that sum to the east/west distance between point O and point xy
			double ry; // r is the point on the same station as o (60) and the same line as xy xy, r, o form a right triangle

			// if the user has specified +lon_0 or +k0 for some reason,
			// we're going to ignore it so that xy is consistent with point O
			lp.lam=lp.lam+lam0;
			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<=EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			xy.x=lp.lam;
			xy.y=-Math.Log(Proj.pj_tsfn(lp.phi, Math.Sin(lp.phi), e)); // Mercator transform xy
			
			oy=-Math.Log(Proj.pj_tsfn(PT_O_PHI, Math.Sin(PT_O_PHI), e));
			
			l1=(xy.y-oy)*Math.Tan(ROTATION_ANGLE);
			l2=-xy.x-l1+PT_O_LAMBDA;
			
			ry=l2*Math.Cos(ROTATION_ANGLE)*Math.Sin(ROTATION_ANGLE)+xy.y;
			ry=Proj.pj_phi2(ctx, Math.Exp(-ry), e); // inverse Mercator
		
			xy.x=PT_O_LINE-Proj.RAD_TO_DEG*(ry-PT_O_PHI)*DEG_TO_LINE/Math.Cos(ROTATION_ANGLE);
			xy.y=PT_O_STATION+Proj.RAD_TO_DEG*(ry-lp.phi)*DEG_TO_STATION/Math.Sin(ROTATION_ANGLE);

			// set a = 1, x0 = 0, and y0 = 0 so that no further unit adjustments are done
			a=1;
			x0=0;
			y0=0;

			return xy;
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;
	
			lp.lam=lp.lam+lam0;
			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<=EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			xy.x=lp.lam;
			xy.y=Math.Log(Math.Tan(Proj.FORTPI+0.5*lp.phi));
			
			double oy=Math.Log(Math.Tan(Proj.FORTPI+0.5*PT_O_PHI));
			
			double l1=(xy.y-oy)*Math.Tan(ROTATION_ANGLE);
			double l2=-xy.x-l1+PT_O_LAMBDA;
			
			double ry=l2*Math.Cos(ROTATION_ANGLE)*Math.Sin(ROTATION_ANGLE)+xy.y;
			ry=Proj.HALFPI-2.0*Math.Atan(Math.Exp(-ry));
			
			xy.x=PT_O_LINE-Proj.RAD_TO_DEG*(ry-PT_O_PHI)*DEG_TO_LINE/Math.Cos(ROTATION_ANGLE);
			xy.y=PT_O_STATION+Proj.RAD_TO_DEG*(ry-lp.phi)*DEG_TO_STATION/Math.Sin(ROTATION_ANGLE);

			a=1;
			x0=0;
			y0=0;

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double ry; // y value of point r
			double oymctr; // Mercator-transformed y value of point O
			double rymctr; // Mercator-transformed ry
			double xymctr; // Mercator-transformed xy.y

			// turn x and y back into Line/Station
			xy.x/=ra;
			xy.y/=ra;

			ry=PT_O_PHI-LINE_TO_RAD*(xy.x-PT_O_LINE)*Math.Cos(ROTATION_ANGLE);

			lp.phi=ry-STATION_TO_RAD*(xy.y-PT_O_STATION)*Math.Sin(ROTATION_ANGLE);

			oymctr=-Math.Log(Proj.pj_tsfn(PT_O_PHI, Math.Sin(PT_O_PHI), e));
			rymctr=-Math.Log(Proj.pj_tsfn(ry, Math.Sin(ry), e));
			xymctr=-Math.Log(Proj.pj_tsfn(lp.phi, Math.Sin(lp.phi), e));

			double l1=(xymctr-oymctr)*Math.Tan(ROTATION_ANGLE);
			double l2=(rymctr-xymctr)/(Math.Cos(ROTATION_ANGLE)*Math.Sin(ROTATION_ANGLE));
			
			lp.lam=PT_O_LAMBDA-(l1+l2);
			
			over=true;
			
			return lp;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.x/=ra;
			xy.y/=ra;

			double ry=PT_O_PHI-LINE_TO_RAD*(xy.x-PT_O_LINE)*Math.Cos(ROTATION_ANGLE);
	
			lp.phi=ry-STATION_TO_RAD*(xy.y-PT_O_STATION)*Math.Sin(ROTATION_ANGLE);

			double oymctr=Math.Log(Math.Tan(Proj.FORTPI+.5*PT_O_PHI));
			double rymctr=Math.Log(Math.Tan(Proj.FORTPI+.5*ry));
			double xymctr=Math.Log(Math.Tan(Proj.FORTPI+.5*lp.phi));

			double l1=(xymctr-oymctr)*Math.Tan(ROTATION_ANGLE);
			double l2=(rymctr-xymctr)/(Math.Cos(ROTATION_ANGLE)*Math.Sin(ROTATION_ANGLE));

			lp.lam=PT_O_LAMBDA-(l1+l2);

			over=true;

			return lp;
		}

		public override PJ Init()
		{
			if(es!=0)
			{ // ellipsoid
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{ // sphere
				inv=s_inverse;
				fwd=s_forward;
			}

			return this;
		}
	}
}

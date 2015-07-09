using System;

namespace Free.Ports.Proj4.Projections
{
	// The Natural Earth projection was designed by Tom Patterson, US National Park
	// Service, in 2007, using Flex Projector. The shape of the original projection
	// was defined at every 5 degrees and piece-wise cubic spline interpolation was
	// used to compute the complete graticule.
	// The code here uses polynomial functions instead of cubic splines and
	// is therefore much simpler to program. The polynomial approximation was
	// developed by Bojan Savric, in collaboration with Tom Patterson and Bernhard
	// Jenny, Institute of Cartography, ETH Zurich. It slightly deviates from
	// Patterson's original projection by adding additional curvature to meridians
	// where they meet the horizontal pole line. This improvement is by intention
	// and designed in collaboration with Tom Patterson.
	// Port to PROJ.4 by Bernhard Jenny, 6 June 2011
	class PJ_natearth : PJ
	{
		public override string Name { get { return "natearth"; } }
		public override string DescriptionName { get { return "Natural Earth"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double A0=0.8707;
		const double A1=-0.131979;
		const double A2=-0.013791;
		const double A3=0.003971;
		const double A4=-0.001529;
		const double B0=1.007226;
		const double B1=0.015085;
		const double B2=-0.044475;
		const double B3=0.028874;
		const double B4=-0.005916;
		const double C0=B0;
		const double C1=(3*B1);
		const double C2=(7*B2);
		const double C3=(9*B3);
		const double C4=(11*B4);
		const double EPS=1e-11;
		const double MAX_Y=(0.8707*0.52*Math.PI);

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double phi2=lp.phi*lp.phi;
			double phi4=phi2*phi2;
			xy.x=lp.lam*(A0+phi2*(A1+phi2*(A2+phi4*phi2*(A3+phi2*A4))));
			xy.y=lp.phi*(B0+phi2*(B1+phi4*(B2+B3*phi2+B4*phi4)));
			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			// make sure y is inside valid range
			if(xy.y>MAX_Y)
			{
				xy.y=MAX_Y;
			}
			else if(xy.y<-MAX_Y)
			{
				xy.y=-MAX_Y;
			}

			// latitude
			double yc=xy.y;
			for(; ; )
			{ // Newton-Raphson
				double y2=yc*yc;
				double y4=y2*y2;
				double f=(yc*(B0+y2*(B1+y4*(B2+B3*y2+B4*y4))))-xy.y;
				double fder=C0+y2*(C1+y4*(C2+C3*y2+C4*y4));
				double tol=f/fder;
				yc-=tol;
				if(Math.Abs(tol)<EPS) break;
			}
			lp.phi=yc;

			// longitude
			{
				double y2=yc*yc;
				lp.lam=xy.x/(A0+y2*(A1+y2*(A2+y2*y2*y2*(A3+y2*A4))));
			}

			return lp;
		}

		public override PJ Init()
		{
			es=0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}

using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_gnom : PJ
	{
		protected double sinph0, cosph0;
		protected gnom_mode mode;

		public override string Name { get { return "gnom"; } }
		public override string DescriptionName { get { return "Gnomonic"; } }
		public override string DescriptionType { get { return "Azi, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double EPS10=1.0e-10;

		protected enum gnom_mode
		{
			N_POLE=0,
			S_POLE=1,
			EQUIT=2,
			OBLIQ=3
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sinphi=Math.Sin(lp.phi);
			double cosphi=Math.Cos(lp.phi);
			double coslam=Math.Cos(lp.lam);

			switch(mode)
			{
				case gnom_mode.EQUIT: xy.y=cosphi*coslam; break;
				case gnom_mode.OBLIQ: xy.y=sinph0*sinphi+cosph0*cosphi*coslam; break;
				case gnom_mode.S_POLE: xy.y=-sinphi; break;
				case gnom_mode.N_POLE: xy.y=sinphi; break;
			}

			if(xy.y<=EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			xy.x=(xy.y=1.0/xy.y)*cosphi*Math.Sin(lp.lam);

			switch(mode)
			{
				case gnom_mode.EQUIT: xy.y*=sinphi; break;
				case gnom_mode.OBLIQ: xy.y*=cosph0*sinphi-sinph0*cosphi*coslam; break;
				case gnom_mode.N_POLE: coslam=-coslam; goto case gnom_mode.S_POLE;
				case gnom_mode.S_POLE: xy.y*=cosphi*coslam; break;
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double rh=Libc.hypot(xy.x, xy.y);
			lp.phi=Math.Atan(rh);
			double sinz=Math.Sin(lp.phi);
			double cosz=Math.Sqrt(1.0-sinz*sinz);

			if(Math.Abs(rh)<=EPS10)
			{
				lp.phi=phi0;
				lp.lam=0.0;
			}
			else
			{
				switch(mode)
				{
					case gnom_mode.OBLIQ:
						lp.phi=cosz*sinph0+xy.y*sinz*cosph0/rh;
						if(Math.Abs(lp.phi)>=1.0) lp.phi=lp.phi>0.0?Proj.HALFPI:-Proj.HALFPI;
						else lp.phi=Math.Asin(lp.phi);
						xy.y=(cosz-sinph0*Math.Sin(lp.phi))*rh;
						xy.x*=sinz*cosph0;
						break;
					case gnom_mode.EQUIT:
						lp.phi=xy.y*sinz/rh;
						if(Math.Abs(lp.phi)>=1.0) lp.phi=lp.phi>0.0?Proj.HALFPI:-Proj.HALFPI;
						else lp.phi=Math.Asin(lp.phi);
						xy.y=cosz*rh;
						xy.x*=sinz;
						break;
					case gnom_mode.S_POLE: lp.phi-=Proj.HALFPI; break;
					case gnom_mode.N_POLE: lp.phi=Proj.HALFPI-lp.phi; xy.y=-xy.y; break;
				}

				lp.lam=Math.Atan2(xy.x, xy.y);
			}

			return lp;
		}

		public override PJ Init()
		{
			if(Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<EPS10) mode=phi0<0.0?gnom_mode.S_POLE:gnom_mode.N_POLE;
			else if(Math.Abs(phi0)<EPS10) mode=gnom_mode.EQUIT;
			else
			{
				mode=gnom_mode.OBLIQ;
				sinph0=Math.Sin(phi0);
				cosph0=Math.Cos(phi0);
			}

			inv=s_inverse;
			fwd=s_forward;
			es=0.0;

			return this;
		}
	}
}

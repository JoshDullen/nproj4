using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_ortho : PJ
	{
		protected double sinph0, cosph0;
		protected ortho_mode mode;

		public override string Name { get { return "ortho"; } }
		public override string DescriptionName { get { return "Orthographic"; } }
		public override string DescriptionType { get { return "Azi, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double EPS10=1.0e-10;

		protected enum ortho_mode
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

			double cosphi=Math.Cos(lp.phi);
			double coslam=Math.Cos(lp.lam);

			switch(mode)
			{
				case ortho_mode.EQUIT:
					if(cosphi*coslam<-EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					xy.y=Math.Sin(lp.phi);
					break;
				case ortho_mode.OBLIQ:
					double sinphi=Math.Sin(lp.phi);
					if(sinph0*sinphi+cosph0*cosphi*coslam<-EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					xy.y=cosph0*sinphi-sinph0*cosphi*coslam;
					break;
				case ortho_mode.N_POLE:
					coslam=-coslam;
					if(Math.Abs(lp.phi-phi0)-EPS10>Proj.HALFPI) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					xy.y=cosphi*coslam;
					break;
				case ortho_mode.S_POLE:
					if(Math.Abs(lp.phi-phi0)-EPS10>Proj.HALFPI) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					xy.y=cosphi*coslam;
					break;
			}
			xy.x=cosphi*Math.Sin(lp.lam);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double rh=Libc.hypot(xy.x, xy.y);
			double sinc=rh;

			if(sinc>1.0)
			{
				if((sinc-1.0)>EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				sinc=1.0;
			}

			double cosc=Math.Sqrt(1.0-sinc*sinc); // in this range OK
			if(Math.Abs(rh)<=EPS10)
			{
				lp.phi=phi0;
				lp.lam=0.0;
			}
			else
			{
				switch(mode)
				{
					case ortho_mode.N_POLE: xy.y=-xy.y; lp.phi=Math.Acos(sinc); break;
					case ortho_mode.S_POLE: lp.phi=-Math.Acos(sinc); break;
					case ortho_mode.EQUIT:
						lp.phi=xy.y*sinc/rh;
						xy.x*=sinc;
						xy.y=cosc*rh;
						if(Math.Abs(lp.phi)>=1.0) lp.phi=lp.phi<0.0?-Proj.HALFPI:Proj.HALFPI;
						else lp.phi=Math.Asin(lp.phi);
						break;
					case ortho_mode.OBLIQ:
						lp.phi=cosc*sinph0+xy.y*sinc*cosph0/rh;
						xy.y=(cosc-sinph0*lp.phi)*rh;
						xy.x*=sinc*cosph0;
						if(Math.Abs(lp.phi)>=1.0) lp.phi=lp.phi<0.0?-Proj.HALFPI:Proj.HALFPI;
						else lp.phi=Math.Asin(lp.phi);
						break;
				}
				lp.lam=(xy.y==0.0&&(mode==ortho_mode.OBLIQ||mode==ortho_mode.EQUIT))?
					(xy.x==0.0?0.0:xy.x<0.0?-Proj.HALFPI:Proj.HALFPI):Math.Atan2(xy.x, xy.y);
			}

			return lp;
		}

		public override PJ Init()
		{
			if(Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<=EPS10) mode=phi0<0.0?ortho_mode.S_POLE:ortho_mode.N_POLE;
			else if(Math.Abs(phi0)>EPS10)
			{
				mode=ortho_mode.OBLIQ;
				sinph0=Math.Sin(phi0);
				cosph0=Math.Cos(phi0);
			}
			else mode=ortho_mode.EQUIT;

			inv=s_inverse;
			fwd=s_forward;
			es=0.0;

			return this;
		}
	}
}

using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_stere : PJ
	{
		protected double phits, sinX1, cosX1, akm1;
		protected stere_mode mode;

		public override string Name { get { return "stere"; } }
		public override string DescriptionName { get { return "Stereographic"; } }
		public override string DescriptionType { get { return "Azi, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "lat_ts="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				if(phits!=Proj.HALFPI) ret.AppendFormat(nc, " +lat_ts={0}", phits*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		const double EPS10=1.0e-10;
		const double TOL=1.0e-8;
		const int NITER=8;
		const double CONV=1e-10;

		protected enum stere_mode
		{
			S_POLE=0,
			N_POLE=1,
			OBLIQ=2,
			EQUIT=3
		}

		static double ssfn_(double phit, double sinphi, double eccen)
		{
			sinphi*=eccen;
			return Math.Tan(0.5*(Proj.HALFPI+phit))*Math.Pow((1.0-sinphi)/(1.0+sinphi), 0.5*eccen);
		}

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sinX=0.0, cosX=0.0;

			double coslam=Math.Cos(lp.lam);
			double sinlam=Math.Sin(lp.lam);
			double sinphi=Math.Sin(lp.phi);

			if(mode==stere_mode.OBLIQ||mode==stere_mode.EQUIT)
			{
				double X=2.0*Math.Atan(ssfn_(lp.phi, sinphi, e))-Proj.HALFPI;
				sinX=Math.Sin(X);
				cosX=Math.Cos(X);
			}

			switch(mode)
			{
				case stere_mode.OBLIQ:
					{
						double A=akm1/(cosX1*(1.0+sinX1*sinX+cosX1*cosX*coslam));
						xy.y=A*(cosX1*sinX-sinX1*cosX*coslam);
						xy.x=A*cosX;
					}
					break;
				case stere_mode.EQUIT:
					{
						double A=2.0*akm1/(1.0+cosX*coslam);
						xy.y=A*sinX;
						xy.x=A*cosX;
					}
					break;
				case stere_mode.S_POLE:
					lp.phi=-lp.phi;
					coslam=-coslam;
					sinphi=-sinphi;
					goto case stere_mode.N_POLE;
				case stere_mode.N_POLE:
					xy.x=akm1*Proj.pj_tsfn(lp.phi, sinphi, e);
					xy.y=-xy.x*coslam;
					break;
			}

			xy.x=xy.x*sinlam;

			return xy;
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sinphi=Math.Sin(lp.phi);
			double cosphi=Math.Cos(lp.phi);
			double coslam=Math.Cos(lp.lam);
			double sinlam=Math.Sin(lp.lam);

			switch(mode)
			{
				case stere_mode.EQUIT:
					xy.y=1.0+cosphi*coslam;
					if(xy.y<=EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					xy.y=akm1/xy.y;
					xy.x=xy.y*cosphi*sinlam;
					xy.y*=sinphi;
					break;
				case stere_mode.OBLIQ:
					xy.y=1.0+sinX1*sinphi+cosX1*cosphi*coslam;
					if(xy.y<=EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					xy.y=akm1/xy.y;
					xy.x=xy.y*cosphi*sinlam;
					xy.y*=cosX1*sinphi-sinX1*cosphi*coslam;
					break;
				case stere_mode.N_POLE:
					coslam=-coslam;
					lp.phi=-lp.phi;
					goto case stere_mode.S_POLE;
				case stere_mode.S_POLE:
					if(Math.Abs(lp.phi-Proj.HALFPI)<TOL) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					xy.y=akm1*Math.Tan(Proj.FORTPI+0.5*lp.phi);
					xy.x=sinlam*xy.y;
					xy.y*=coslam;
					break;
			}

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double sinphi, tp=0.0, phi_l=0.0, halfe=0.0, halfpi=0.0;

			double rho=Libc.hypot(xy.x, xy.y);
			switch(mode)
			{
				case stere_mode.OBLIQ:
				case stere_mode.EQUIT:
					tp=2.0*Math.Atan2(rho*cosX1, akm1);
					double cosphi=Math.Cos(tp);
					sinphi=Math.Sin(tp);
					if(rho==0.0) phi_l=Math.Asin(cosphi*sinX1);
					else phi_l=Math.Asin(cosphi*sinX1+(xy.y*sinphi*cosX1/rho));
					tp=Math.Tan(0.5*(Proj.HALFPI+phi_l));
					xy.x*=sinphi;
					xy.y=rho*cosX1*cosphi-xy.y*sinX1*sinphi;
					halfpi=Proj.HALFPI;
					halfe=0.5*e;
					break;
				case stere_mode.N_POLE:
					xy.y=-xy.y;
					goto case stere_mode.S_POLE;
				case stere_mode.S_POLE:
					tp=-rho/akm1;
					phi_l=Proj.HALFPI-2.0*Math.Atan(tp);
					halfpi=-Proj.HALFPI;
					halfe=-0.5*e;
					break;
			}

			for(int i=NITER-1; i>=0; i++)
			{
				sinphi=e*Math.Sin(phi_l);
				lp.phi=2.0*Math.Atan(tp*Math.Pow((1.0+sinphi)/(1.0-sinphi), halfe))-halfpi;
				if(Math.Abs(phi_l-lp.phi)<CONV)
				{
					if(mode==stere_mode.S_POLE) lp.phi=-lp.phi;
					lp.lam=(xy.x==0.0&&xy.y==0.0)?0.0:Math.Atan2(xy.x, xy.y);
					return lp;
				}
				phi_l=lp.phi;
			}

			Proj.pj_ctx_set_errno(ctx, -20);
			return lp;

		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double rh=Libc.hypot(xy.x, xy.y);
			double c=2.0*Math.Atan(rh/akm1);
			double sinc=Math.Sin(c);
			double cosc=Math.Cos(c);

			lp.lam=0.0;
			switch(mode)
			{
				case stere_mode.EQUIT:
					if(Math.Abs(rh)<=EPS10) lp.phi=0.0;
					else lp.phi=Math.Asin(xy.y*sinc/rh);
					if(cosc!=0.0||xy.x!=0.0) lp.lam=Math.Atan2(xy.x*sinc, cosc*rh);
					break;
				case stere_mode.OBLIQ:
					if(Math.Abs(rh)<=EPS10) lp.phi=phi0;
					else lp.phi=Math.Asin(cosc*sinX1+xy.y*sinc*cosX1/rh);
					c=cosc-sinX1*Math.Sin(lp.phi);
					if(c!=0.0||xy.x!=0.0) lp.lam=Math.Atan2(xy.x*sinc*cosX1, c*rh);
					break;
				case stere_mode.N_POLE:
					xy.y=-xy.y;
					goto case stere_mode.S_POLE;
				case stere_mode.S_POLE:
					if(Math.Abs(rh)<=EPS10) lp.phi=phi0;
					else lp.phi=Math.Asin(mode==stere_mode.S_POLE?-cosc:cosc);
					lp.lam=(xy.x==0.0&&xy.y==0.0)?0.0:Math.Atan2(xy.x, xy.y);
					break;
			}

			return lp;
		}

		protected PJ_stere setup()
		{ // general initialization
			double t=Math.Abs(phi0);

			if(Math.Abs(t-Proj.HALFPI)<EPS10) mode=phi0<0.0?stere_mode.S_POLE:stere_mode.N_POLE;
			else mode=t>EPS10?stere_mode.OBLIQ:stere_mode.EQUIT;

			phits=Math.Abs(phits);
			if(es!=0)
			{
				switch(mode)
				{
					case stere_mode.N_POLE:
					case stere_mode.S_POLE:
						if(Math.Abs(phits-Proj.HALFPI)<EPS10) akm1=2.0*k0/Math.Sqrt(Math.Pow(1+e, 1+e)*Math.Pow(1-e, 1-e));
						else
						{
							t=Math.Sin(phits);
							akm1=Math.Cos(phits)/Proj.pj_tsfn(phits, t, e);
							t*=e;
							akm1/=Math.Sqrt(1.0-t*t);
						}
						break;
					case stere_mode.EQUIT:
					case stere_mode.OBLIQ:
						t=Math.Sin(phi0);
						double X=2.0*Math.Atan(ssfn_(phi0, t, e))-Proj.HALFPI;
						t*=e;
						akm1=2.0*k0*Math.Cos(phi0)/Math.Sqrt(1.0-t*t);
						sinX1=Math.Sin(X);
						cosX1=Math.Cos(X);
						break;
				}
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{
				switch(mode)
				{
					case stere_mode.OBLIQ:
						sinX1=Math.Sin(phi0);
						cosX1=Math.Cos(phi0);
						goto case stere_mode.EQUIT;
					case stere_mode.EQUIT:
						akm1=2.0*k0;
						break;
					case stere_mode.S_POLE:
					case stere_mode.N_POLE:
						akm1=Math.Abs(phits-Proj.HALFPI)>=EPS10?Math.Cos(phits)/Math.Tan(Proj.FORTPI-0.5*phits):2.0*k0;
						break;
				}
				inv=s_inverse;
				fwd=s_forward;
			}

			return this;
		}

		public override PJ Init()
		{
			phits=Proj.pj_param_t(ctx, parameters, "lat_ts")?Proj.pj_param_r(ctx, parameters, "lat_ts"):Proj.HALFPI;

			return setup();
		}
	}

	class PJ_ups : PJ_stere
	{
		public override string Name { get { return "ups"; } }
		public override string DescriptionName { get { return "Universal Polar Stereographic"; } }
		public override string DescriptionParameters { get { return "south"; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				if(phi0<0) ret.Append(" +south");
				return ret.ToString();
			}
		}

		public override PJ Init()
		{
			// International Ellipsoid
			phi0=Proj.pj_param_b(ctx, parameters, "south")?-Proj.HALFPI:Proj.HALFPI;
			if(es==0) { Proj.pj_ctx_set_errno(ctx, -34); return null; }
			k0=0.994;
			x0=2000000.0;
			y0=2000000.0;
			phits=Proj.HALFPI;
			lam0=0.0;

			return setup();
		}
	}
}

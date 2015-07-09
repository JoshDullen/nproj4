using System;
using System.Text;
using Free.Ports.Proj4.FactorsDeriv;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_lcc : PJ
	{
		protected double phi1, phi2, n, rho0, c;
		protected bool ellips;

		public override string Name { get { return "lcc"; } }
		public override string DescriptionName { get { return "Lambert Conformal Conic"; } }
		public override string DescriptionType { get { return "Conic, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "lat_1= and lat_2= or lat_0"; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", phi1*Proj.RAD_TO_DEG);
				if(phi2==phi1) ret.AppendFormat(nc, " +lat_2={0}", phi2*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		public void_LP_PJ_FACTORS spc { get; protected set; }

		const double EPS10=1.0e-10;

		// ellipsoid & spheroid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;
			double rho=0.0;

			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<EPS10)
			{
				if((lp.phi*n)<=0.0) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
				rho=0.0;
			}
			else rho=c*(ellips?Math.Pow(Proj.pj_tsfn(lp.phi, Math.Sin(lp.phi), e), n):Math.Pow(Math.Tan(Proj.FORTPI+0.5*lp.phi), -n));

			lp.lam*=n;
			xy.x=k0*(rho*Math.Sin(lp.lam));
			xy.y=k0*(rho0-rho*Math.Cos(lp.lam));

			return xy;
		}

		// ellipsoid & spheroid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;
			double rho=0.0;

			xy.x/=k0;
			xy.y/=k0;
			xy.y=rho0-xy.y;
			rho=Libc.hypot(xy.x, xy.y);
			if(rho!=0.0)
			{
				if(n<0.0)
				{
					rho=-rho;
					xy.x=-xy.x;
					xy.y=-xy.y;
				}
				if(ellips)
				{
					lp.phi=Proj.pj_phi2(ctx, Math.Pow(rho/c, 1.0/n), e);
					if(lp.phi==Libc.HUGE_VAL) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				}
				else lp.phi=2.0*Math.Atan(Math.Pow(c/rho, 1.0/n))-Proj.HALFPI;
				lp.lam=Math.Atan2(xy.x, xy.y)/n;
			}
			else
			{
				lp.lam=0.0;
				lp.phi=n>0.0?Proj.HALFPI:-Proj.HALFPI;
			}

			return lp;
		}

		void fac(LP lp, ref FACTORS fac)
		{
			double rho=0.0;

			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<EPS10)
			{
				if((lp.phi*n)<=0.0) return;
				rho=0.0;
			}
			else rho=c*(ellips?Math.Pow(Proj.pj_tsfn(lp.phi, Math.Sin(lp.phi), e), n):Math.Pow(Math.Tan(Proj.FORTPI+0.5*lp.phi), -n));
			fac.code|=Factors.IS_ANAL_HK+Factors.IS_ANAL_CONV;
			fac.k=fac.h=k0*n*rho/Proj.pj_msfn(Math.Sin(lp.phi), Math.Cos(lp.phi), es);
			fac.conv=-n*lp.lam;
		}

		public override PJ Init()
		{
			phi1=Proj.pj_param_r(ctx, parameters, "lat_1");
			if(Proj.pj_param_t(ctx, parameters, "lat_2")) phi2=Proj.pj_param_r(ctx, parameters, "lat_2");
			else
			{
				phi2=phi1;
				if(!Proj.pj_param_t(ctx, parameters, "lat_0")) phi0=phi1;
			}

			if(Math.Abs(phi1+phi2)<EPS10) { Proj.pj_ctx_set_errno(ctx, -21); return null; }
			double sinphi=n=Math.Sin(phi1);
			double cosphi=Math.Cos(phi1);
			bool secant=Math.Abs(phi1-phi2)>=EPS10;

			ellips=(es!=0.0);
			if(ellips)
			{
				e=Math.Sqrt(es);
				double m1=Proj.pj_msfn(sinphi, cosphi, es);
				double ml1=Proj.pj_tsfn(phi1, sinphi, e);
				if(secant)
				{ // secant cone
					sinphi=Math.Sin(phi2);
					n=Math.Log(m1/Proj.pj_msfn(sinphi, Math.Cos(phi2), es));
					n/=Math.Log(ml1/Proj.pj_tsfn(phi2, sinphi, e));
				}
				c=rho0=m1*Math.Pow(ml1, -n)/n;
				rho0*=(Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<EPS10)?0.0:Math.Pow(Proj.pj_tsfn(phi0, Math.Sin(phi0), e), n);
			}
			else
			{
				if(secant) n=Math.Log(cosphi/Math.Cos(phi2))/Math.Log(Math.Tan(Proj.FORTPI+0.5*phi2)/Math.Tan(Proj.FORTPI+0.5*phi1));
				c=cosphi*Math.Pow(Math.Tan(Proj.FORTPI+0.5*phi1), n)/n;
				rho0=(Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<EPS10)?0.0:c*Math.Pow(Math.Tan(Proj.FORTPI+0.5*phi0), -n);
			}

			inv=e_inverse;
			fwd=e_forward;
			spc=fac;

			return this;
		}
	}
}

using System;
using Free.Ports.Proj4.LibcStuff;

// projection scale factors

namespace Free.Ports.Proj4.FactorsDeriv
{
	public static partial class Factors
	{
		public const int IS_ANAL_XL_YL=01;	// derivatives of lon analytic
		public const int IS_ANAL_XP_YP=02;	// derivatives of lat analytic
		public const int IS_ANAL_HK=04;		// h and k analytic
		public const int IS_ANAL_CONV=010;	// convergence analytic

		public static bool pj_factors(LP lp, PJ P, double h, ref FACTORS fac)
		{
			IFactors PFac=P as IFactors;

			const double EPS=1.0e-12;
			const double DEFAULT_H=1e-5;

			// check for forward and latitude or longitude overange
			double t=Math.Abs(lp.phi)-Proj.HALFPI;
			if(t>EPS||Math.Abs(lp.lam)>10.0)
			{
				Proj.pj_ctx_set_errno(P.ctx, -14);
				return true;
			}

			// proceed
			Libc.errno=Proj.pj_errno=0;
			P.ctx.last_errno=0;
			if(h<EPS) h=DEFAULT_H;
			if(Math.Abs(lp.phi)>(Proj.HALFPI-h))
				lp.phi=lp.phi<0.0?(-Proj.HALFPI+h):(Proj.HALFPI-h); // adjust to value around pi/2 where derived still exists
			else if(P.geoc) lp.phi=Math.Atan(P.rone_es*Math.Tan(lp.phi));

			lp.lam-=P.lam0; // compute del lp.lam
			if(!P.over) lp.lam=Proj.adjlon(lp.lam); // adjust del longitude
			if(PFac!=null) PFac.spc(lp, ref fac); // get what projection analytic values

			if((fac.code&(IS_ANAL_XL_YL+IS_ANAL_XP_YP))!=(IS_ANAL_XL_YL+IS_ANAL_XP_YP))
			{
				DERIVS der;
				if(pj_deriv(lp, h, P, out der)) return true;

				if((fac.code&IS_ANAL_XL_YL)==0)
				{
					fac.der.x_l=der.x_l;
					fac.der.y_l=der.y_l;
				}
				if((fac.code&IS_ANAL_XP_YP)==0)
				{
					fac.der.x_p=der.x_p;
					fac.der.y_p=der.y_p;
				}
			}

			double cosphi=Math.Cos(lp.phi);

			double r;
			if((fac.code&IS_ANAL_HK)==0)
			{
				fac.h=Libc.hypot(fac.der.x_p, fac.der.y_p);
				fac.k=Libc.hypot(fac.der.x_l, fac.der.y_l)/cosphi;
				if(P.es!=0)
				{
					t=Math.Sin(lp.phi);
					t=1.0-P.es*t*t;
					double n=Math.Sqrt(t);
					fac.h*=t*n/P.one_es;
					fac.k*=n;
					r=t*t/P.one_es;
				}
				else r=1.0;
			}
			else if(P.es!=0)
			{
				r=Math.Sin(lp.phi);
				r=1.0-P.es*r*r;
				r=r*r/P.one_es;
			}
			else r=1.0;

			// convergence
			if((fac.code&IS_ANAL_CONV)==0)
			{
				fac.conv=-Math.Atan2(fac.der.y_l, fac.der.x_l);
				if((fac.code&IS_ANAL_XL_YL)!=0) fac.code|=IS_ANAL_CONV;
			}

			// areal scale factor
			fac.s=(fac.der.y_p*fac.der.x_l-fac.der.x_p*fac.der.y_l)*r/cosphi;

			// meridian-parallel angle theta prime
			fac.thetap=Proj.aasin(P.ctx, fac.s/(fac.h*fac.k));

			// Tissot ellips axis
			t=fac.k*fac.k+fac.h*fac.h;
			fac.a=Math.Sqrt(t+2.0*fac.s);
			t=(t=t-2.0*fac.s)<=0.0?0.0:Math.Sqrt(t);
			fac.b=0.5*(fac.a-t);
			fac.a=0.5*(fac.a+t);

			// omega
			fac.omega=2.0*Proj.aasin(P.ctx, (fac.a-fac.b)/(fac.a+fac.b));

			return false;
		}
	}
}

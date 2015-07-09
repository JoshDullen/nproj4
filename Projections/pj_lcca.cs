// PROJ.4 Cartographic Projection System

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_lcca : PJ
	{
		protected double r0, l, M0, C;
		protected double[] en;

		public override string Name { get { return "lcca"; } }
		public override string DescriptionName { get { return "Lambert Conformal Conic Alternative"; } }
		public override string DescriptionType { get { return "Conic, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "lat_0="; } }
		public override bool Invertible { get { return true; } }

		const int MAX_ITER=10;
		const double DEL_TOL=1e-12;

		// func to compute dr
		static double fS(double S, double C) { return S*(1.0+S*S*C); }

		// deriv of fs
		static double fSp(double S, double C) { return 1.0+3.0*S*S*C; }

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double S=Proj.pj_mlfn(lp.phi, Math.Sin(lp.phi), Math.Cos(lp.phi), en)-M0;
			double dr=fS(S, C);
			double r=r0-dr;

			lp.lam*=l;
			xy.x=k0*r*Math.Sin(lp.lam);
			xy.y=k0*(r0-r*Math.Cos(lp.lam));

			return xy;
		}

		// ellipsoid & spheroid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.x/=k0;
			xy.y/=k0;

			double theta=Math.Atan2(xy.x, r0-xy.y);
			double dr=xy.y-xy.x*Math.Tan(0.5*theta);
			lp.lam=theta/l;
			double S=dr;

			int i=MAX_ITER;
			for(; i>0; i--)
			{
				double dif=(fS(S, C)-dr)/fSp(S, C);
				S-=dif;
				if(Math.Abs(dif)<DEL_TOL) break;
			}
			if(i==0) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

			lp.phi=Proj.pj_inv_mlfn(ctx, S+M0, es, en);

			return lp;
		}

		public override PJ Init()
		{
			en=Proj.pj_enfn(es);
			if(en==null) return null;

			if(!Proj.pj_param_t(ctx, parameters, "lat_0")) { Proj.pj_ctx_set_errno(ctx, 50); return null; }
			if(phi0==0.0) { Proj.pj_ctx_set_errno(ctx, 51); return null; }

			l=Math.Sin(phi0);
			M0=Proj.pj_mlfn(phi0, l, Math.Cos(phi0), en);
			double s2p0=l*l;
			double R0=1.0/(1.0-es*s2p0);
			double N0=Math.Sqrt(R0);
			R0*=one_es*N0;
			double tan0=Math.Tan(phi0);
			r0=N0/tan0;
			C=1.0/(6.0*R0*N0);

			inv=e_inverse;
			fwd=e_forward;

			return this;
		}
	}
}

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_somerc : PJ
	{
		protected double K, c, hlf_e, kR, cosp0, sinp0;

		public override string Name { get { return "somerc"; } }
		public override string DescriptionName { get { return "Swiss. Obl. Mercator (CH1903)"; } }
		public override string DescriptionType { get { return "Cyl, Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double EPS=1e-10;
		const int NITER=6;

		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sp=e*Math.Sin(lp.phi);
			double phip=2.0*Math.Atan(Math.Exp(c*(Math.Log(Math.Tan(Proj.FORTPI+0.5*lp.phi))-hlf_e*Math.Log((1.0+sp)/(1.0-sp)))+K))-Proj.HALFPI;
			double lamp=c*lp.lam;
			double cp=Math.Cos(phip);
			double phipp=Proj.aasin(ctx, cosp0*Math.Sin(phip)-sinp0*cp*Math.Cos(lamp));
			double lampp=Proj.aasin(ctx, cp*Math.Sin(lamp)/Math.Cos(phipp));
			xy.x=kR*lampp;
			xy.y=kR*Math.Log(Math.Tan(Proj.FORTPI+0.5*phipp));

			return xy;
		}

		// ellipsoid & spheroid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double phipp=2.0*(Math.Atan(Math.Exp(xy.y/kR))-Proj.FORTPI);
			double lampp=xy.x/kR;
			double cp=Math.Cos(phipp);
			double phip=Proj.aasin(ctx, cosp0*Math.Sin(phipp)+sinp0*cp*Math.Cos(lampp));
			double lamp=Proj.aasin(ctx, cp*Math.Sin(lampp)/Math.Cos(phip));
			double con=(K-Math.Log(Math.Tan(Proj.FORTPI+0.5*phip)))/c;

			int i=NITER;
			for(; i>0; i--)
			{
				double esp=e*Math.Sin(phip);
				double delp=(con+Math.Log(Math.Tan(Proj.FORTPI+0.5*phip))-hlf_e*Math.Log((1.0+esp)/(1.0-esp)))*(1.0-esp*esp)*Math.Cos(phip)*rone_es;
				phip-=delp;
				if(Math.Abs(delp)<EPS) break;
			}

			if(i==0) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
			lp.phi=phip;
			lp.lam=lamp/c;

			return lp;
		}

		public override PJ Init()
		{
			hlf_e=0.5*e;
			double cp=Math.Cos(phi0);
			cp*=cp;
			c=Math.Sqrt(1+es*cp*cp*rone_es);
			double sp=Math.Sin(phi0);
			sinp0=sp/c;
			double phip0=Proj.aasin(ctx, sinp0);
			cosp0=Math.Cos(phip0);
			sp*=e;
			K=Math.Log(Math.Tan(Proj.FORTPI+0.5*phip0))-c*(Math.Log(Math.Tan(Proj.FORTPI+0.5*phi0))-hlf_e*Math.Log((1.0+sp)/(1.0-sp)));
			kR=k0*Math.Sqrt(one_es)/(1.0-sp*sp);
			inv=e_inverse;
			fwd=e_forward;

			return this;
		}
	}
}

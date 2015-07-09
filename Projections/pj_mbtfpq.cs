using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_mbtfpq : PJ
	{
		public override string Name { get { return "mbtfpq"; } }
		public override string DescriptionName { get { return "McBryde-Thomas Flat-Polar Quartic"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const int NITER=20;
		const double EPS=1.0e-7;
		const double ONETOL=1.0000001;
		const double C=1.70710678118654752440;
		const double RC=0.58578643762690495119;
		const double FYC=1.87475828462269495505;
		const double RYC=0.53340209679417701685;
		const double FXC=0.31245971410378249250;
		const double RXC=3.20041258076506210122;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double c=C*Math.Sin(lp.phi);
			for(int i=NITER; i>0; i--)
			{
				double th1=(Math.Sin(.5*lp.phi)+Math.Sin(lp.phi)-c)/(0.5*Math.Cos(0.5*lp.phi)+Math.Cos(lp.phi));
				lp.phi-=th1;
				if(Math.Abs(th1)<EPS) break;
			}

			xy.x=FXC*lp.lam*(1.0+2.0*Math.Cos(lp.phi)/Math.Cos(0.5*lp.phi));
			xy.y=FYC*Math.Sin(0.5*lp.phi);
			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double t;

			lp.phi=RYC*xy.y;
			if(Math.Abs(lp.phi)>1.0)
			{
				if(Math.Abs(lp.phi)>ONETOL) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				else if(lp.phi<0.0) { t=-1.0; lp.phi=-Proj.PI; }
				else { t=1.0; lp.phi=Proj.PI; }
			}
			else
			{
				t=lp.phi;
				lp.phi=2.0*Math.Asin(t);
			}

			lp.lam=RXC*xy.x/(1.0+2.0*Math.Cos(lp.phi)/Math.Cos(0.5*lp.phi));
			lp.phi=RC*(t+Math.Sin(lp.phi));
			if(Math.Abs(lp.phi)>1.0)
			{
				if(Math.Abs(lp.phi)>ONETOL) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				else lp.phi=lp.phi<0.0?-Proj.HALFPI:Proj.HALFPI;
			}
			else lp.phi=Math.Asin(lp.phi);

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

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_nell_h : PJ
	{
		public override string Name { get { return "nell_h"; } }
		public override string DescriptionName { get { return "Nell-Hammer"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const int NITER=9;
		const double EPS=1.0e-7;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=0.5*lp.lam*(1.0+Math.Cos(lp.phi));
			xy.y=2.0*(lp.phi-Math.Tan(0.5*lp.phi));

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			int i;
			double p=0.5*xy.y;

			for(i=NITER; i>0; i--)
			{
				double c=Math.Cos(0.5*lp.phi);
				double V=(lp.phi-Math.Tan(lp.phi/2)-p)/(1.0-0.5/(c*c));
				lp.phi-=V;
				if(Math.Abs(V)<EPS) break;
			}

			if(i==0)
			{
				lp.phi=p<0.0?-Proj.HALFPI:Proj.HALFPI;
				lp.lam=2.0*xy.x;
			}
			else lp.lam=2.0*xy.x/(1.0+Math.Cos(lp.phi));

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

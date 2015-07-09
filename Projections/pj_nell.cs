using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_nell : PJ
	{
		public override string Name { get { return "nell"; } }
		public override string DescriptionName { get { return "Nell"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const int MAX_ITER=10;
		const double LOOP_TOL=1.0e-7;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double k=2.0*Math.Sin(lp.phi);
			double V=lp.phi*lp.phi;
			lp.phi*=1.00371+V*(-0.0935382+V*-0.011412);

			for(int i=MAX_ITER; i>0; i--)
			{
				V=(lp.phi+Math.Sin(lp.phi)-k)/(1.0+Math.Cos(lp.phi));
				lp.phi-=V;
				if(Math.Abs(V)<LOOP_TOL) break;
			}

			xy.x=0.5*lp.lam*(1.0+Math.Cos(lp.phi));
			xy.y=lp.phi;
			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.lam=2.0*xy.x/(1.0+Math.Cos(xy.y));
			lp.phi=Proj.aasin(ctx, 0.5*(xy.y+Math.Sin(xy.y)));

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

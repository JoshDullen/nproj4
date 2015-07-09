using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_fouc_s : PJ
	{
		protected double n, n1;

		public override string Name { get { return "fouc_s"; } }
		public override string DescriptionName { get { return "Foucaut Sinusoidal"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return "n="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +n={0}", n);
				return ret.ToString();
			}
		}

		const int MAX_ITER=10;
		const double LOOP_TOL=1e-7;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double t=Math.Cos(lp.phi);
			xy.x=lp.lam*t/(n+n1*t);
			xy.y=n*lp.phi+n1*Math.Sin(lp.phi);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			if(n!=0)
			{
				lp.phi=xy.y;
				int i=MAX_ITER;
				for(; i>0; i--)
				{
					double v=(n*lp.phi+n1*Math.Sin(lp.phi)-xy.y)/(n+n1*Math.Cos(lp.phi));
					lp.phi-=v;
					if(Math.Abs(v)<LOOP_TOL) break;
				}

				if(i==0) lp.phi=xy.y<0.0?-Proj.HALFPI:Proj.HALFPI;
			}
			else lp.phi=Proj.aasin(ctx, xy.y);

			double V=Math.Cos(lp.phi);
			lp.lam=xy.x*(n+n1*V)/V;

			return lp;
		}

		public override PJ Init()
		{
			n=Proj.pj_param_d(ctx, parameters, "n");
			if(n<0.0||n>1.0) { Proj.pj_ctx_set_errno(ctx, -99); return null; }
			n1=1.0-n;
			es=0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}

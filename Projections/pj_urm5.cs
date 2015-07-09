using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_urm5 : PJ
	{
		protected double m, rmn, q3, n;

		public override string Name { get { return "urm5"; } }
		public override string DescriptionName { get { return "Urmaev V"; } }
		public override string DescriptionType { get { return "PCyl, Sph, no inv."; } }
		public override string DescriptionParameters { get { return "n= q= [m= or alpha=]"; } }
		public override bool Invertible { get { return false; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +n={0}", n);
				ret.AppendFormat(nc, " +q={0}", q3*3);
				ret.AppendFormat(nc, " +m={0}", m);
				return ret.ToString();
			}
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double t=Proj.aasin(ctx, n*Math.Sin(lp.phi));
			lp.phi=t;
			xy.x=m*lp.lam*Math.Cos(lp.phi);
			t*=t;
			xy.y=lp.phi*(1.0+t*q3)*rmn;

			return xy;
		}

		public override PJ Init()
		{
			n=Proj.pj_param_d(ctx, parameters, "n");
			q3=Proj.pj_param_d(ctx, parameters, "q")/3.0;
			if(Proj.pj_param_t(ctx, parameters, "m")&&!Proj.pj_param_t(ctx, parameters, "alpha")) m=Proj.pj_param_d(ctx, parameters, "m");
			else
			{
				double alpha=Proj.pj_param_r(ctx, parameters, "alpha");
				double t=n*Math.Sin(alpha);
				m=Math.Cos(alpha)/Math.Sqrt(1.0-t*t);
			}
			rmn=1.0/(m*n);
			es=0.0;
			inv=null;
			fwd=s_forward;

			return this;
		}
	}
}

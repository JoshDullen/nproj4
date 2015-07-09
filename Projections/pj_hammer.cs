using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_hammer : PJ
	{
		protected double w, m, rm;

		public override string Name { get { return "hammer"; } }
		public override string DescriptionName { get { return "Hammer & Eckert-Greifendorff"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return "W= M="; } }
		public override bool Invertible { get { return false; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				if(w!=0.5) ret.AppendFormat(nc, " +W={0}", w);
				if(rm!=1.0) ret.AppendFormat(nc, " +M={0}", 1/rm);
				return ret.ToString();
			}
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double cosphi=Math.Cos(lp.phi);
			lp.lam*=w;
			double d=Math.Sqrt(2.0/(1.0+cosphi*Math.Cos(lp.lam)));
			xy.x=m*d*cosphi*Math.Sin(lp.lam);
			xy.y=rm*d*Math.Sin(lp.phi);

			return xy;
		}

		public override PJ Init()
		{
			if(Proj.pj_param_t(ctx, parameters, "W"))
			{
				w=Math.Abs(Proj.pj_param_d(ctx, parameters, "W"));
				if(w<=0.0) { Proj.pj_ctx_set_errno(ctx, -27); return null; }
			}
			else w=0.5;

			if(Proj.pj_param_t(ctx, parameters, "M"))
			{
				m=Math.Abs(Proj.pj_param_d(ctx, parameters, "M"));
				if(m<=0.0) { Proj.pj_ctx_set_errno(ctx, -27); return null; }
			}
			else m=1.0;

			rm=1.0/m;
			m/=w;
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}
}

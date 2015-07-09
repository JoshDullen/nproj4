using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_wink1 : PJ
	{
		protected double cosphi1;

		public override string Name { get { return "wink1"; } }
		public override string DescriptionName { get { return "Winkel I"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return "lat_ts="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_ts={0}", Math.Acos(cosphi1)*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=0.5*lp.lam*(cosphi1+Math.Cos(lp.phi));
			xy.y=lp.phi;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y;
			lp.lam=2.0*xy.x/(cosphi1+Math.Cos(lp.phi));

			return lp;
		}

		public override PJ Init()
		{
			cosphi1=Math.Cos(Proj.pj_param_r(ctx, parameters, "lat_ts"));
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}

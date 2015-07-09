using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_eqc : PJ
	{
		protected double rc;

		public override string Name { get { return "eqc"; } }
		public override string DescriptionName { get { return "Equidistant Cylindrical (Plate Caree)"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return "lat_ts= [lat_0=0]"; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_ts={0}", Math.Acos(rc)*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=rc*lp.lam;
			xy.y=lp.phi-phi0;
			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.lam=xy.x/rc;
			lp.phi=xy.y+phi0;
			return lp;
		}

		public override PJ Init()
		{
			rc=Math.Cos(Proj.pj_param_r(ctx, parameters, "lat_ts"));
			if(rc<=0.0) { Proj.pj_ctx_set_errno(ctx, -24); return null; }
			inv=s_inverse;
			fwd=s_forward;
			es=0.0;

			return this;
		}
	}
}

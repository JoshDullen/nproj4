using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_tcc : PJ
	{
		public override string Name { get { return "tcc"; } }
		public override string DescriptionName { get { return "Transverse Central Cylindrical"; } }
		public override string DescriptionType { get { return "Cyl, Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const double EPS10=1.0e-10;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double b=Math.Cos(lp.phi)*Math.Sin(lp.lam);
			double bt=1.0-b*b;
			if(bt<EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			xy.x=b/Math.Sqrt(bt);
			xy.y=Math.Atan2(Math.Tan(lp.phi), Math.Cos(lp.lam));

			return xy;
		}

		public override PJ Init()
		{
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}
}

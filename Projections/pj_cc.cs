using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_cc : PJ
	{
		public override string Name { get { return "cc"; } }
		public override string DescriptionName { get { return "Central Cylindrical"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double EPS10=1.0e-10;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<=EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
			xy.x=lp.lam;
			xy.y=Math.Tan(lp.phi);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=Math.Atan(xy.y);
			lp.lam=xy.x;

			return lp;
		}

		public override PJ Init()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}

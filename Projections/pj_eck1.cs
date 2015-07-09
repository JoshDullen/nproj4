using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_eck1 : PJ
	{
		public override string Name { get { return "eck1"; } }
		public override string DescriptionName { get { return "Eckert I"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double FC=0.92131773192356127802;
		const double RP=0.31830988618379067154;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=FC*lp.lam*(1.0-RP*Math.Abs(lp.phi));
			xy.y=FC*lp.phi;
			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y/FC;
			lp.lam=xy.x/(FC*(1.0-RP*Math.Abs(lp.phi)));
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

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_larr : PJ
	{
		public override string Name { get { return "larr"; } }
		public override string DescriptionName { get { return "Larrivee"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const double SIXTH=0.16666666666666666;

		// sphere
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=0.5*lp.lam*(1.0+Math.Sqrt(Math.Cos(lp.phi)));
			xy.y=lp.phi/(Math.Cos(0.5*lp.phi)*Math.Cos(SIXTH*lp.lam));

			return xy;
		}

		public override PJ Init()
		{
			fwd=s_forward;
			inv=null;
			es=0.0;

			return this;
		}
	}
}

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_eck5 : PJ
	{
		public override string Name { get { return "eck5"; } }
		public override string DescriptionName { get { return "Eckert V"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double XF=0.44101277172455148219;
		const double RXF=2.26750802723822639137;
		const double YF=0.88202554344910296438;
		const double RYF=1.13375401361911319568;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=XF*(1.0+Math.Cos(lp.phi))*lp.lam;
			xy.y=YF*lp.phi;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=RYF*xy.y;
			lp.lam=RXF*xy.x/(1.0+Math.Cos(lp.phi));

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

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_fahey : PJ
	{
		public override string Name { get { return "fahey"; } }
		public override string DescriptionName { get { return "Fahey"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double TOL=1.0e-6;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=Math.Tan(0.5*lp.phi);
			xy.y=1.819152*xy.x;
			xy.x=0.819152*lp.lam*Proj.asqrt(1-xy.x*xy.x);
			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y/=1.819152;
			lp.phi=2.0*Math.Atan(xy.y);
			xy.y=1.0-xy.y*xy.y;
			lp.lam=Math.Abs(xy.y)<TOL?0.0:xy.x/(0.819152*Math.Sqrt(xy.y));

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

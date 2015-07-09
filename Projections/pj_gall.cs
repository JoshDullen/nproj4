using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_gall : PJ
	{
		public override string Name { get { return "gall"; } }
		public override string DescriptionName { get { return "Gall (Gall Stereographic)"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double YF=1.70710678118654752440;
		const double XF=0.70710678118654752440;
		const double RYF=0.58578643762690495119;
		const double RXF=1.41421356237309504880;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=XF*lp.lam;
			xy.y=YF*Math.Tan(0.5*lp.phi);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.lam=RXF*xy.x;
			lp.phi=2.0*Math.Atan(xy.y*RYF);

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

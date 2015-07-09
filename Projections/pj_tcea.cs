using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_tcea : PJ
	{
		protected double rk0;

		public override string Name { get { return "tcea"; } }
		public override string DescriptionName { get { return "Transverse Cylindrical Equal Area"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=rk0*Math.Cos(lp.phi)*Math.Sin(lp.lam);
			xy.y=k0*(Math.Atan2(Math.Tan(lp.phi), Math.Cos(lp.lam))-phi0);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y=xy.y*rk0+phi0;
			xy.x*=k0;
			double t=Math.Sqrt(1.0-xy.x*xy.x);
			lp.phi=Math.Asin(t*Math.Sin(xy.y));
			lp.lam=Math.Atan2(xy.x, t*Math.Cos(xy.y));

			return lp;
		}

		public override PJ Init()
		{
			rk0=1/k0;
			inv=s_inverse;
			fwd=s_forward;
			es=0.0;

			return this;
		}
	}
}

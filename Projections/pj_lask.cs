namespace Free.Ports.Proj4.Projections
{
	class PJ_lask : PJ
	{
		public override string Name { get { return "lask"; } }
		public override string DescriptionName { get { return "Laskowski"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const double a10=0.975534;
		const double a12=-0.119161;
		const double a32=-0.0143059;
		const double a14=-0.0547009;
		const double b01=1.00384;
		const double b21=0.0802894;
		const double b03=0.0998909;
		const double b41=0.000199025;
		const double b23=-0.0285500;
		const double b05=-0.0491032;

		// sphere
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double l2=lp.lam*lp.lam;
			double p2=lp.phi*lp.phi;
			xy.x=lp.lam*(a10+p2*(a12+l2*a32+p2*a14));
			xy.y=lp.phi*(b01+l2*(b21+p2*b23+l2*b41)+p2*(b03+p2*b05));

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

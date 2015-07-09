namespace Free.Ports.Proj4.Projections
{
	class PJ_gins8 : PJ
	{
		public override string Name { get { return "gins8"; } }
		public override string DescriptionName { get { return "Ginsburg VIII (TsNIIGAiK)"; } }
		public override string DescriptionType { get { return "PCyl, Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const double Cl=0.000952426;
		const double Cp=0.162388;
		const double C12=0.08333333333333333;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double t=lp.phi*lp.phi;
			xy.y=lp.phi*(1.0+t*C12);
			xy.x=lp.lam*(1.0-Cp*t);
			t=lp.lam*lp.lam;
			xy.x*=(0.87-Cl*t*t);

			return xy;
		}

		public override PJ Init()
		{
			es=0.0;
			inv=null;
			fwd=s_forward;

			return this;
		}
	}
}

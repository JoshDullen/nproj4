namespace Free.Ports.Proj4.Projections
{
	class PJ_eck3 : PJ
	{
		protected double C_x, C_y, A, B;

		public override string Name { get { return "eck3"; } }
		public override string DescriptionName { get { return "Eckert III"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.y=C_y*lp.phi;
			xy.x=C_x*lp.lam*(A+Proj.asqrt(1.0-B*lp.phi*lp.phi));
			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y/C_y;
			lp.lam=xy.x/(C_x*(A+Proj.asqrt(1.0-B*lp.phi*lp.phi)));
			return lp;
		}

		protected PJ setup()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}

		public override PJ Init()
		{
			C_x=0.42223820031577120149;
			C_y=0.84447640063154240298;
			A=1.0;
			B=0.4052847345693510857755;

			return setup();
		}
	}

	class PJ_kav7 : PJ_eck3
	{
		public override string Name { get { return "kav7"; } }
		public override string DescriptionName { get { return "Kavraisky VII"; } }

		public override PJ Init()
		{
			//C_x=0.2632401569273184856851; ???
			C_x=0.8660254037844;
			C_y=1.0;
			A=0.0;
			B=0.30396355092701331433;

			return setup();
		}
	}

	class PJ_wag6 : PJ_eck3
	{
		public override string Name { get { return "wag6"; } }
		public override string DescriptionName { get { return "Wagner VI"; } }

		public override PJ Init()
		{
			C_x=C_y=0.94745;
			A=0.0;
			B=0.30396355092701331433;

			return setup();
		}
	}

	class PJ_putp1 : PJ_eck3
	{
		public override string Name { get { return "putp1"; } }
		public override string DescriptionName { get { return "Putnins P1"; } }

		public override PJ Init()
		{
			C_x=1.89490;
			C_y=0.94745;
			A=-0.5;
			B=0.30396355092701331433;

			return setup();
		}
	}
}

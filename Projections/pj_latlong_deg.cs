namespace Free.Ports.Proj4.Projections
{
	class PJ_latlong_deg : PJ
	{
		public override string Name { get { return "latlong_deg"; } }
		public override string DescriptionName { get { return "Lat/long in degree (Geodetic)"; } }
		public override string DescriptionType { get { return ""; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		XY forward(LP lp)
		{
			XY xy;
			xy.x=lp.lam*Proj.RAD_TO_DEG/a;
			xy.y=lp.phi*Proj.RAD_TO_DEG/a;
			return xy;
		}

		LP inverse(XY xy)
		{
			LP lp;
			lp.phi=xy.y*Proj.DEG_TO_RAD*a;
			lp.lam=xy.x*Proj.DEG_TO_RAD*a;
			return lp;
		}

		public override PJ Init()
		{
			is_latlong=false;
			x0=0.0;
			y0=0.0;
			inv=inverse;
			fwd=forward;

			return this;
		}
	}
}

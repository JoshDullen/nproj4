using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_gstmerc : PJ
	{
		protected double phic, c, n1, n2, XS, YS;

		public override string Name { get { return "gstmerc"; } }
		public override string DescriptionName { get { return "Gauss-Schreiber Transverse Mercator (aka Gauss-Laborde Reunion)"; } }
		public override string DescriptionType { get { return "Cyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double L=n1*lp.lam;
			double Ls=c+n1*Math.Log(Proj.pj_tsfn(-1.0*lp.phi, -1.0*Math.Sin(lp.phi), e));
			double sinLs1=Math.Sin(L)/Math.Cosh(Ls);
			double Ls1=Math.Log(Proj.pj_tsfn(-1.0*Math.Asin(sinLs1), 0.0, 0.0));
			xy.x=(XS+n2*Ls1)*ra;
			xy.y=(YS+n2*Math.Atan(Math.Sinh(Ls)/Math.Cos(L)))*ra;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double L=Math.Atan(Math.Sinh((xy.x*a-XS)/n2)/Math.Cos((xy.y*a-YS)/n2));
			double sinC=Math.Sin((xy.y*a-YS)/n2)/Math.Cosh((xy.x*a-XS)/n2);
			double LC=Math.Log(Proj.pj_tsfn(-1.0*Math.Asin(sinC), 0.0, 0.0));
			lp.lam=L/n1;
			lp.phi=-1.0*Proj.pj_phi2(ctx, Math.Exp((LC-c)/n1), e);

			return lp;
		}

		public override PJ Init()
		{
			n1=Math.Sqrt(1.0+es*Math.Pow(Math.Cos(phi0), 4.0)/(1.0-es));
			phic=Math.Asin(Math.Sin(phi0)/n1);
			c=Math.Log(Proj.pj_tsfn(-1.0*phic, 0.0, 0.0))-n1*Math.Log(Proj.pj_tsfn(-1.0*phi0, -1.0*Math.Sin(phi0), e));
			n2=k0*a*Math.Sqrt(1.0-es)/(1.0-es*Math.Sin(phi0)*Math.Sin(phi0));
			XS=0; // -P.x0
			YS=-1.0*n2*phic; // -P.y0
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}

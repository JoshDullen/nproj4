using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_wink2 : PJ
	{
		protected double cosphi1;

		public override string Name { get { return "wink2"; } }
		public override string DescriptionName { get { return "Winkel II"; } }
		public override string DescriptionType { get { return "PCyl, Sph, no inv."; } }
		public override string DescriptionParameters { get { return "lat_1="; } }
		public override bool Invertible { get { return false; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", Math.Acos(cosphi1)*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		const double LOOP_TOL=1.0e-7;
		const int MAX_ITER=10;
		const double TWO_D_PI=0.636619772367581343;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.y=lp.phi*TWO_D_PI;
			double k=Proj.PI*Math.Sin(lp.phi);
			lp.phi*=1.8;
			int i=MAX_ITER;
			for(; i>0; i--)
			{
				double V=(lp.phi+Math.Sin(lp.phi)-k)/(1.0+Math.Cos(lp.phi));
				lp.phi-=V;
				if(Math.Abs(V)<LOOP_TOL) break;
			}

			if(i==0) lp.phi=(lp.phi<0.0)?-Proj.HALFPI:Proj.HALFPI;
			else lp.phi*=0.5;

			xy.x=0.5*lp.lam*(Math.Cos(lp.phi)+cosphi1);
			xy.y=Proj.FORTPI*(Math.Sin(lp.phi)+xy.y);

			return xy;
		}

		public override PJ Init()
		{
			cosphi1=Math.Cos(Proj.pj_param_r(ctx, parameters, "lat_1"));
			es=0.0;
			inv=null;
			fwd=s_forward;

			return this;
		}
	}
}

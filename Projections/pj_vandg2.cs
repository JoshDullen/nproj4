using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_vandg2 : PJ
	{
		protected bool vdg3;

		public override string Name { get { return "vandg2"; } }
		public override string DescriptionName { get { return "van der Grinten II"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const double TOL=1e-10;
		const double TWORPI=0.63661977236758134308;

		// spheroid
		protected XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double x1, at;

			double bt=Math.Abs(TWORPI*lp.phi);
			double ct=1.0-bt*bt;
			if(ct<0.0) ct=0.0;
			else ct=Math.Sqrt(ct);

			if(Math.Abs(lp.lam)<TOL)
			{
				xy.x=0.0;
				xy.y=Proj.PI*(lp.phi<0.0?-bt:bt)/(1.0+ct);
			}
			else
			{
				at=0.5*Math.Abs(Proj.PI/lp.lam-lp.lam/Proj.PI);
				if(vdg3)
				{
					x1=bt/(1.0+ct);
					xy.x=Proj.PI*(Math.Sqrt(at*at+1.0-x1*x1)-at);
					xy.y=Proj.PI*x1;
				}
				else
				{
					x1=(ct*Math.Sqrt(1.0+at*at)-at*ct*ct)/(1.0+at*at*bt*bt);
					xy.x=Proj.PI*x1;
					xy.y=Proj.PI*Math.Sqrt(1.0-x1*(x1+2.0*at)+TOL);
				}

				if(lp.lam<0.0) xy.x=-xy.x;
				if(lp.phi<0.0) xy.y=-xy.y;
			}

			return xy;
		}

		public override PJ Init()
		{
			vdg3=false;
			inv=null;
			fwd=s_forward;
			
			return this;
		}
	}

	class PJ_vandg3 : PJ_vandg2
	{
		public override string Name { get { return "vandg3"; } }
		public override string DescriptionName { get { return "van der Grinten III"; } }

		public override PJ Init()
		{
			vdg3=true;
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}
}

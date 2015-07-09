using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_vandg4 : PJ
	{
		public override string Name { get { return "vandg4"; } }
		public override string DescriptionName { get { return "van der Grinten IV"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const double TOL=1e-10;
		const double TWORPI=0.63661977236758134308;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;
			
			if(Math.Abs(lp.phi)<TOL)
			{
				xy.x=lp.lam;
				xy.y=0.0;
			}
			else if(Math.Abs(lp.lam)<TOL||Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<TOL)
			{
				xy.x=0.0;
				xy.y=lp.phi;
			}
			else
			{
				double bt=Math.Abs(TWORPI*lp.phi);
				double bt2=bt*bt;
				double ct=0.5*(bt*(8.0-bt*(2.0+bt2))-5.0)/(bt2*(bt-1.0));
				double ct2=ct*ct;
				double dt=TWORPI*lp.lam;
				dt=dt+1.0/dt;
				dt=Math.Sqrt(dt*dt-4.0);
				if((Math.Abs(lp.lam)-Proj.HALFPI)<0.0) dt=-dt;
				double dt2=dt*dt;
				double x1=bt+ct;
				x1*=x1;
				double t=bt+3.0*ct;
				double ft=x1*(bt2+ct2*dt2-1.0)+(1.0-bt2)*(bt2*(t*t+4.0*ct2)+ct2*(12.0*bt*ct+4.0*ct2));
				x1=(dt*(x1+ct2-1.0)+2.0*Math.Sqrt(ft))/(4.0*x1+dt2);
				xy.x=Proj.HALFPI*x1;
				xy.y=Proj.HALFPI*Math.Sqrt(1.0+dt*Math.Abs(x1)-x1*x1);
				if(lp.lam<0.0) xy.x=-xy.x;
				if(lp.phi<0.0) xy.y=-xy.y;
			}

			return xy;
		}

		public override PJ Init()
		{
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}
}

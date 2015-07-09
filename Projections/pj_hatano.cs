using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_hatano : PJ
	{
		public override string Name { get { return "hatano"; } }
		public override string DescriptionName { get { return "Hatano Asymmetrical Equal Area"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const int NITER=20;
		const double EPS=1.0e-7;
		const double ONETOL=1.000001;
		const double CN=2.67595;
		const double CS=2.43763;
		const double RCN=0.37369906014686373063;
		const double RCS=0.41023453108141924738;
		const double FYCN=1.75859;
		const double FYCS=1.93052;
		const double RYCN=0.56863737426006061674;
		const double RYCS=0.51799515156538134803;
		const double FXC=0.85;
		const double RXC=1.17647058823529411764;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			int i=NITER;

			double c=Math.Sin(lp.phi)*(lp.phi<0.0?CS:CN);
			for(; i>0; i--)
			{
				double th1=(lp.phi+Math.Sin(lp.phi)-c)/(1.0+Math.Cos(lp.phi));
				lp.phi-=th1;
				if(Math.Abs(th1)<EPS) break;
			}

			lp.phi*=0.5;

			xy.x=FXC*lp.lam*Math.Cos(lp.phi);
			xy.y=Math.Sin(lp.phi)*(lp.phi<0.0?FYCS:FYCN);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double th=xy.y*(xy.y<0.0?RYCS:RYCN);

			if(Math.Abs(th)>1.0)
			{
				if(Math.Abs(th)>ONETOL) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				else th=th>0.0?Proj.HALFPI:-Proj.HALFPI;
			}
			else th=Math.Asin(th);

			lp.lam=RXC*xy.x/Math.Cos(th);
			th+=th;
			lp.phi=(th+Math.Sin(th))*(xy.y<0.0?RCS:RCN);

			if(Math.Abs(lp.phi)>1.0)
			{
				if(Math.Abs(lp.phi)>ONETOL) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				else lp.phi=lp.phi>0.0?Proj.HALFPI:-Proj.HALFPI;
			}
			else lp.phi=Math.Asin(lp.phi);

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

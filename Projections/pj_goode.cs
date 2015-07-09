using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_goode : PJ
	{
		protected PJ sinu, moll;

		public override string Name { get { return "goode"; } }
		public override string DescriptionName { get { return "Goode Homolosine"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double Y_COR=0.05280;
		const double PHI_LIM=0.71093078197902358062;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(Math.Abs(lp.phi)<=PHI_LIM) xy=sinu.fwd(lp);
			else
			{
				xy=moll.fwd(lp);
				xy.y-=(lp.phi>=0.0?Y_COR:-Y_COR);
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			if(Math.Abs(xy.y)<=PHI_LIM) lp=sinu.inv(xy);
			else
			{
				xy.y+=(xy.y>=0.0?Y_COR:-Y_COR);
				lp=moll.inv(xy);
			}

			return lp;
		}

		public override PJ Init()
		{
			es=0.0;

			try
			{
				sinu=new PJ_sinu();
				moll=new PJ_moll();
			}
			catch
			{
				return null;
			}

			sinu.es=0;
			sinu.ctx=ctx;
			moll.ctx=ctx;

			sinu=sinu.Init();
			moll=moll.Init();
			if(sinu==null||moll==null) return null;

			fwd=s_forward;
			inv=s_inverse;

			return this;
		}
	}
}

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_poly : PJ
	{
		protected double ml0;
		protected double[] en;

		public override string Name { get { return "poly"; } }
		public override string DescriptionName { get { return "Polyconic (American)"; } }
		public override string DescriptionType { get { return "Conic, Sph&Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double TOL=1.0e-10;
		const double CONV=1.0e-10;
		const int N_ITER=10;
		const int I_ITER=20;
		const double ITOL=1e-12;

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(Math.Abs(lp.phi)<=TOL) { xy.x=lp.lam; xy.y=-ml0; }
			else
			{
				double sp=Math.Sin(lp.phi);
				double cp=Math.Cos(lp.phi);
				double ms=Math.Abs(cp)>TOL?Proj.pj_msfn(sp, cp, es)/sp:0.0;
				lp.lam*=sp;
				xy.x=ms*Math.Sin(lp.lam);
				xy.y=(Proj.pj_mlfn(lp.phi, sp, cp, en)-ml0)+ms*(1.0-Math.Cos(lp.lam));
			}

			return xy;
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(Math.Abs(lp.phi)<=TOL) { xy.x=lp.lam; xy.y=ml0; }
			else
			{
				double cot=1.0/Math.Tan(lp.phi);
				double E=lp.lam*Math.Sin(lp.phi);
				xy.x=Math.Sin(E)*cot;
				xy.y=lp.phi-phi0+cot*(1.0-Math.Cos(E));
			}

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y+=ml0;

			if(Math.Abs(xy.y)<=TOL) { lp.lam=xy.x; lp.phi=0.0; return lp; }

			double c;
			double r=xy.y*xy.y+xy.x*xy.x;
			lp.phi=xy.y;

			int i=I_ITER;
			for(; i>0; i--)
			{
				double sp=Math.Sin(lp.phi);
				double cp=Math.Cos(lp.phi);
				double s2ph=sp*cp;
				if(Math.Abs(cp)<ITOL) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

				double mlp=Math.Sqrt(1.0-es*sp*sp);
				c=sp*mlp/cp;
				double ml=Proj.pj_mlfn(lp.phi, sp, cp, en);
				double mlb=ml*ml+r;
				mlp=one_es/(mlp*mlp*mlp);

				double dPhi=(ml+ml+c*mlb-2.0*xy.y*(c*ml+1.0))/(es*s2ph*(mlb-2.0*xy.y*ml)/c+2.0*(xy.y-ml)*(c*mlp-1.0/s2ph)-mlp-mlp);
				lp.phi+=dPhi;

				if(Math.Abs(dPhi)<=ITOL) break;
			}

			if(i==0) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

			c=Math.Sin(lp.phi);
			lp.lam=Math.Asin(xy.x*Math.Tan(lp.phi)*Math.Sqrt(1.0-es*c*c))/Math.Sin(lp.phi);

			return lp;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y=phi0+xy.y;
			if(Math.Abs(xy.y)<=TOL) { lp.lam=xy.x; lp.phi=0.0; }
			else
			{

				lp.phi=xy.y;
				double B=xy.x*xy.x+xy.y*xy.y;
				int i=N_ITER;
				double dphi;
				do
				{
					double tp=Math.Tan(lp.phi);
					dphi=(xy.y*(lp.phi*tp+1.0)-lp.phi-0.5*(lp.phi*lp.phi+B)*tp)/((lp.phi-xy.y)/tp-1.0);
					lp.phi-=dphi;
					i--;
				} while(Math.Abs(dphi)>CONV&&i!=0);

				if(i==0) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

				lp.lam=Math.Asin(xy.x*Math.Tan(lp.phi))/Math.Sin(lp.phi);
			}

			return lp;
		}

		public override PJ Init()
		{
			if(es!=0)
			{
				en=Proj.pj_enfn(es);
				if(en==null) return null;
				ml0=Proj.pj_mlfn(phi0, Math.Sin(phi0), Math.Cos(phi0), en);
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{
				ml0=-phi0;
				inv=s_inverse;
				fwd=s_forward;
			}

			return this;
		}
	}
}

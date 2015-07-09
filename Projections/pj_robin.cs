using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_robin : PJ
	{
		public override string Name { get { return "robin"; } }
		public override string DescriptionName { get { return "Robinson"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// note: following terms based upon 5 deg. intervals in degrees.
		//
		// Some background on these coefficients is available at:
		// http://article.gmane.org/gmane.comp.gis.proj-4.devel/6039
		// http://trac.osgeo.org/proj/ticket/113
		struct COEFS
		{
			public double c0, c1, c2, c3;

			public COEFS(double c0, double c1, double c2, double c3)
			{
				this.c0=c0;
				this.c1=c1;
				this.c2=c2;
				this.c3=c3;
			}

			public double V(double z)
			{
				return c0+z*(c1+z*(c2+z*c3));
			}

			public double DV(double z)
			{
				return c1+z*(c2+c2+z*3.0*c3);
			}
		}

		static COEFS[] X=new COEFS[]
		{
			new COEFS(1.0000, 2.2199e-17, -7.15515e-05, 3.1103e-06),
			new COEFS(0.9986, -0.000482243, -2.4897e-05, -1.3309e-06),
			new COEFS(0.9954, -0.00083103, -4.48605e-05, -9.86701e-07),
			new COEFS(0.9900, -0.00135364, -5.9661e-05, 3.6777e-06),
			new COEFS(0.9822, -0.00167442, -4.49547e-06, -5.72411e-06),
			new COEFS(0.9730, -0.00214868, -9.03571e-05, 1.8736e-08),
			new COEFS(0.9600, -0.00305085, -9.00761e-05, 1.64917e-06),
			new COEFS(0.9427, -0.00382792, -6.53386e-05, -2.6154e-06),
			new COEFS(0.9216, -0.00467746, -0.00010457, 4.81243e-06),
			new COEFS(0.8962, -0.00536223, -3.23831e-05, -5.43432e-06),
			new COEFS(0.8679, -0.00609363, -0.000113898, 3.32484e-06),
			new COEFS(0.8350, -0.00698325, -6.40253e-05, 9.34959e-07),
			new COEFS(0.7986, -0.00755338, -5.00009e-05, 9.35324e-07),
			new COEFS(0.7597, -0.00798324, -3.5971e-05, -2.27626e-06),
			new COEFS(0.7186, -0.00851367, -7.01149e-05, -8.6303e-06),
			new COEFS(0.6732, -0.00986209, -0.000199569, 1.91974e-05),
			new COEFS(0.6213, -0.010418, 8.83923e-05, 6.24051e-06),
			new COEFS(0.5722, -0.00906601, 0.000182, 6.24051e-06),
			new COEFS(0.5322, -0.00677797, 0.000275608, 6.24051e-06),
		};

		static COEFS[] Y=new COEFS[]
		{
			new COEFS(-5.20417e-18, 0.0124, 1.21431e-18, -8.45284e-11),
			new COEFS(0.0620, 0.0124, -1.26793e-09, 4.22642e-10),
			new COEFS(0.1240, 0.0124, 5.07171e-09, -1.60604e-09),
			new COEFS(0.1860, 0.0123999, -1.90189e-08, 6.00152e-09),
			new COEFS(0.2480, 0.0124002, 7.10039e-08, -2.24e-08),
			new COEFS(0.3100, 0.0123992, -2.64997e-07, 8.35986e-08),
			new COEFS(0.3720, 0.0124029, 9.88983e-07, -3.11994e-07),
			new COEFS(0.4340, 0.0123893, -3.69093e-06, -4.35621e-07),
			new COEFS(0.4958, 0.0123198, -1.02252e-05, -3.45523e-07),
			new COEFS(0.5571, 0.0121916, -1.54081e-05, -5.82288e-07),
			new COEFS(0.6176, 0.0119938, -2.41424e-05, -5.25327e-07),
			new COEFS(0.6769, 0.011713, -3.20223e-05, -5.16405e-07),
			new COEFS(0.7346, 0.0113541, -3.97684e-05, -6.09052e-07),
			new COEFS(0.7903, 0.0109107, -4.89042e-05, -1.04739e-06),
			new COEFS(0.8435, 0.0103431, -6.4615e-05, -1.40374e-09),
			new COEFS(0.8936, 0.00969686, -6.4636e-05, -8.547e-06),
			new COEFS(0.9394, 0.00840947, -0.000192841, -4.2106e-06),
			new COEFS(0.9761, 0.00616527, -0.000256, -4.2106e-06),
			new COEFS(1.0000, 0.00328947, -0.000319159, -4.2106e-06),
		};

		const double FXC=0.8487;
		const double FYC=1.3523;
		const double C1=11.45915590261646417544;
		const double RC1=0.08726646259971647884;
		const int NODES=18;
		const double ONEEPS=1.000001;
		const double EPS=1.0e-8;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double dphi=Math.Abs(lp.phi);
			int i=(int)Math.Floor(dphi*C1);
			if(i>=NODES) i=NODES-1;

			dphi=Proj.RAD_TO_DEG*(dphi-RC1*i);
			xy.x=X[i].V(dphi)*FXC*lp.lam;
			xy.y=Y[i].V(dphi)*FYC;
			if(lp.phi<0.0) xy.y=-xy.y;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.lam=xy.x/FXC;
			lp.phi=Math.Abs(xy.y/FYC);
			if(lp.phi>=1.0)
			{ // simple pathologic cases
				if(lp.phi>ONEEPS) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				else
				{
					lp.phi=xy.y<0.0?-Proj.HALFPI:Proj.HALFPI;
					lp.lam/=X[NODES].c0;
				}
			}
			else
			{ // general problem
				// in Y space, reduce to table interval
				int i=(int)Math.Floor(lp.phi*NODES);

				for(; ; )
				{
					if(Y[i].c0>lp.phi) i--;
					else if(Y[i+1].c0<=lp.phi) i++;
					else break;
				}

				COEFS T=Y[i];

				// first guess, linear interp
				double t=5.0*(lp.phi-T.c0)/(Y[i+1].c0-T.c0);

				// make into root
				T.c0-=lp.phi;
				for(; ; )
				{ // Newton-Raphson reduction
					double t1=T.V(t)/T.DV(t);
					t-=t1;
					if(Math.Abs(t1)<EPS) break;
				}

				lp.phi=(5*i+t)*Proj.DEG_TO_RAD;
				if(xy.y<0.0) lp.phi=-lp.phi;
				lp.lam/=X[i].V(t);
			}

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

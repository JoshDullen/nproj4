using System;
using Free.Ports.Proj4.LibcStuff;

// procedures for evaluating Tseries
namespace Free.Ports.Proj4.Approximation
{
	public static partial class Approximation
	{
		const double NEAR_ONE=1.00001;

		static double ceval(PW_COEF[] C, int n, projUV w, projUV w2)
		{
			double d=0, dd=0;
			int j;
			PW_COEF C1;

			for(; n!=0; n--)
			{
				C1=C[n];
				j=C1.m;
				double tmp=d;
				if(j!=0)
				{
					double vd=0.0, vdd=0.0;
					j--;
					for(; j!=0; j--)
					{
						double tmp2=vd;
						vd=w2.v*tmp2-vdd+C1.c[j];
						vdd=tmp2;
					}
					d=w2.u*tmp-dd+w.v*vd-vdd+0.5*C1.c[j];
				}
				else d=w2.u*tmp-dd;
				dd=tmp;
			}

			C1=C[0];
			j=C1.m;
			if(j!=0)
			{
				double vd=0.0, vdd=0.0;
				j--;
				for(; j!=0; j--)
				{
					double tmp2=vd;
					vd=w2.v*tmp2-vdd+C1.c[j];
					vdd=tmp2;
				}
				return w.u*d-dd+0.5*(w.v*vd-vdd+0.5*C1.c[j]);
			}

			return w.u*d-dd;
		}

		// bivariate Chebyshev polynomial entry point
		public static projUV bcheval(projUV @in, Tseries T)
		{
			projUV w2, w, @out;

			// scale to +-1
			w.u=(@in.u+@in.u-T.a.u)*T.b.u;
			w.v=(@in.v+@in.v-T.a.v)*T.b.v;

			if(Math.Abs(w.u)>NEAR_ONE||Math.Abs(w.v)>NEAR_ONE)
			{
				@out.u=@out.v=Libc.HUGE_VAL;
				Proj.pj_errno=-36;
			}
			else
			{ // double evaluation
				w2.u=w.u+w.u;
				w2.v=w.v+w.v;
				@out.u=ceval(T.cu, T.mu, w, w2);
				@out.v=ceval(T.cv, T.mv, w, w2);
			}

			return @out;
		}

		// bivariate power polynomial entry point
		public static projUV bpseval(projUV @in, Tseries T)
		{
			projUV @out;
			@out.u=@out.v=0.0;

			for(int i=T.mu; i>=0; --i)
			{
				double row=0.0;
				int m=T.cu[i].m;
				if(m!=0)
				{
					int c=m;
					while((m--)!=0) row=T.cu[i].c[--c]+@in.v*row;
				}
				@out.u=row+@in.u*@out.u;
			}

			for(int i=T.mv; i>=0; --i)
			{
				double row=0.0;
				int m=T.cv[i].m;
				if(m!=0)
				{
					int c=m;
					while((m--)!=0) row=T.cv[i].c[--c]+@in.v*row;
				}
				@out.v=row+@in.u*@out.v;
			}

			return @out;
		}

		// general entry point selecting evaluation mode
		public static projUV biveval(projUV @in, Tseries T)
		{
			if(T.power!=0) return bpseval(@in, T);
			return bcheval(@in, T);
		}
	}
}

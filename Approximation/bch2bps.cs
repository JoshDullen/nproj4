// convert bivariate w Chebyshev series to w Power series
namespace Free.Ports.Proj4.Approximation
{
	public static partial class Approximation
	{
		// basic support procedures

		// clear matrix rows to zero
		static void bclear(projUV[][] p, int n, int m)
		{
			for(int i=0; i<n; i++)
				for(int j=0; j<m; j++)
					p[i][j].u=p[i][j].v=0;
		}

		// move vector
		static void bmove(projUV[] a, projUV[] b, int n)
		{
			for(int i=0; i<n; i++) a[i]=b[i];
		}

		// a = m * b - c
		static void submop(projUV[] a, double m, projUV[] b, projUV[] c, int n)
		{
			for(int i=0; i<n; i++)
			{
				a[i].u=m*b[i].u-c[i].u;
				a[i].v=m*b[i].v-c[i].v;
			}
		}

		// a = b - c
		static void subop(projUV[] a, projUV[] b, projUV[] c, int n)
		{
			for(int i=0; i<n; i++)
			{
				a[i].u=b[i].u-c[i].u;
				a[i].v=b[i].v-c[i].v;
			}
		}

		// multiply vector a by scalar m
		static void dmult(projUV[] a, double m, int n)
		{
			for(int i=0; i<n; i++)
			{
				a[i].u*=m;
				a[i].v*=m;
			}
		}

		// row adjust a = a - m * b
		static void dadd(projUV[] a, projUV[] b, double m, int n)
		{
			for(int i=0; i<n; i++)
			{
				a[i].u-=m*b[i].u;
				a[i].v-=m*b[i].v;
			}
		}

		// convert row to pover series
		static void rows(projUV[] c, projUV[] d, int n)
		{
			projUV[] dd=new projUV[n];
			projUV sv;
			sv.u=sv.v=0.0;

			for(int j=0; j<n; j++) d[j]=dd[j]=sv;

			d[0]=c[n-1];
			for(int j=n-2; j>=1; j--)
			{
				for(int k=n-j; k>=1; k--)
				{
					sv=d[k];
					d[k].u=2.0*d[k-1].u-dd[k].u;
					d[k].v=2.0*d[k-1].v-dd[k].v;
					dd[k]=sv;
				}
				sv=d[0];
				d[0].u=-dd[0].u+c[j].u;
				d[0].v=-dd[0].v+c[j].v;
				dd[0]=sv;
			}

			for(int j=n-1; j>=1; j--)
			{
				d[j].u=d[j-1].u-dd[j].u;
				d[j].v=d[j-1].v-dd[j].v;
			}

			d[0].u=-dd[0].u+0.5*c[0].u;
			d[0].v=-dd[0].v+0.5*c[0].v;
		}

		// convert columns to power series
		static void cols(projUV[][] c, projUV[][] d, int nu, int nv)
		{
			projUV[][] dd=new projUV[nu][];
			for(int i=0; i<nu; i++) dd[i]=new projUV[nv];
			projUV[] sv=new projUV[nv];

			bclear(d, nu, nv);
			bclear(dd, nu, nv);
			bmove(d[0], c[nu-1], nv);

			for(int j=nu-2; j>=1; j--)
			{
				for(int k=nu-j; k>=1; k--)
				{
					bmove(sv, d[k], nv);
					submop(d[k], 2.0, d[k-1], dd[k], nv);
					bmove(dd[k], sv, nv);
				}
				bmove(sv, d[0], nv);
				subop(d[0], c[j], dd[0], nv);
				bmove(dd[0], sv, nv);
			}

			for(int j=nu-1; j>=1; j--) subop(d[j], d[j-1], dd[j], nv);
			submop(d[0], 0.5, c[0], dd[0], nv);
		}

		// row adjust for range -1 to 1 to a to b
		static void rowshft(double a, double b, projUV[] d, int n)
		{
			double cnst=2.0/(b-a);
			double fac=cnst;
			for(int j=1; j<n; j++)
			{
				d[j].u*=fac;
				d[j].v*=fac;
				fac*=cnst;
			}

			cnst=0.5*(a+b);
			for(int j=0; j<=n-2; ++j)
			{
				for(int k=n-2; k>=j; k--)
				{
					d[k].u-=cnst*d[k+1].u;
					d[k].v-=cnst*d[k+1].v;
				}
			}
		}

		// column adjust for range -1 to 1 to a to b
		static void colshft(double a, double b, projUV[][] d, int n, int m)
		{
			double cnst=2.0/(b-a);
			double fac=cnst;
			for(int j=1; j<n; j++)
			{
				dmult(d[j], fac, m);
				fac*=cnst;
			}

			cnst=0.5*(a+b);
			for(int j=0; j<=n-2; j++)
				for(int k=n-2; k>=j; k--)
					dadd(d[k], d[k+1], cnst, m);
		}

		// entry point
		public static bool bch2bps(projUV a, projUV b, projUV[][] c, int nu, int nv)
		{
			if(nu<1||nv<1) return false;

			try
			{
				projUV[][] d=new projUV[nu][];
				for(int i=0; i<nu; i++)
				{
					d[i]=new projUV[nv];

					// do rows to power series
					rows(c[i], d[i], nv);
					rowshft(a.v, b.v, d[i], nv);
				}

				// do columns to power series
				cols(d, c, nu, nv);
				colshft(a.u, b.u, c, nu, nv);

				return true;
			}
			catch
			{
				return false;
			}
		}
	}
}

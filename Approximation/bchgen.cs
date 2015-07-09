using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Approximation
{
	public static partial class Approximation
	{
		public delegate projUV ProjFunc(projUV val);

		// generate double bivariate Chebychev polynomial
		public static bool bchgen(projUV a, projUV b, int nu, int nv, projUV[][] f, ProjFunc func)
		{
			projUV arg, bma, bpa;
			projUV[] t, c;

			bma.u=0.5*(b.u-a.u); bma.v=0.5*(b.v-a.v);
			bpa.u=0.5*(b.u+a.u); bpa.v=0.5*(b.v+a.v);
			for(int i=0; i<nu; i++)
			{
				arg.u=Math.Cos(Proj.PI*(i+0.5)/nu)*bma.u+bpa.u;
				for(int j=0; j<nv; j++)
				{
					arg.v=Math.Cos(Proj.PI*(j+0.5)/nv)*bma.v+bpa.v;
					f[i][j]=func(arg);
					if(f[i][j].u==Libc.HUGE_VAL) return true;
				}
			}

			try
			{
				c=new projUV[nu];
			}
			catch
			{
				return true;
			}

			double fac=2.0/nu;
			for(int j=0; j<nv; j++)
			{
				for(int i=0; i<nu; i++)
				{
					arg.u=arg.v=0.0;
					for(int k=0; k<nu; k++)
					{
						double d=Math.Cos(Proj.PI*i*(k+0.5)/nu);
						arg.u+=f[k][j].u*d;
						arg.v+=f[k][j].v*d;
					}
					arg.u*=fac;
					arg.v*=fac;
					c[i]=arg;
				}
				for(int i=0; i<nu; i++) f[i][j]=c[i];
			}

			try
			{
				c=new projUV[nv];
			}
			catch
			{
				return true;
			}

			fac=2.0/nv;
			for(int i=0; i<nu; i++)
			{
				t=f[i];
				for(int j=0; j<nv; j++)
				{
					arg.u=arg.v=0.0;
					for(int k=0; k<nv; k++)
					{
						double d=Math.Cos(Proj.PI*j*(k+0.5)/nv);
						arg.u+=t[k].u*d;
						arg.v+=t[k].v*d;
					}
					arg.u*=fac;
					arg.v*=fac;
					c[j]=arg;
				}
				f[i]=c;
				c=t;
			}

			return false;
		}
	}
}

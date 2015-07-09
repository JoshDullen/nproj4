// This code was written in C by Nathan Wagner
// and ported to C# by the Authors

using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_isea : PJ
	{
		protected isea_dgg dgg;

		public override string Name { get { return "isea"; } }
		public override string DescriptionName { get { return "Icosahedral Snyder Equal Area"; } }
		public override string DescriptionType { get { return "Sph"; } }
		public override string DescriptionParameters { get { return "orient= azi= lon_0= lat_0= aperture= resolution= mode= rescale"; } }
		public override bool Invertible { get { return false; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();

				if(dgg.o_lat==Math.PI/2.0&&dgg.o_lon==0.0) ret.Append(" +orient=pole");

				if(dgg.o_az!=0) ret.AppendFormat(nc, " +azi={0}", dgg.o_az*Proj.RAD_TO_DEG);

				ret.AppendFormat(nc, " +aperture={0}", dgg.aperture);
				ret.AppendFormat(nc, " +resolution={0}", dgg.resolution);

				switch(dgg.output)
				{
					case isea_address_form.ISEA_PLANE: ret.Append(" +mode=plane"); break;
					case isea_address_form.ISEA_Q2DI: ret.Append(" +mode=di"); break;
					case isea_address_form.ISEA_Q2DD: ret.Append(" +mode=dd"); break;
					case isea_address_form.ISEA_HEX: ret.Append(" +mode=hex"); break;
					default: throw new InvalidOperationException("Project settings invalid, or incomplete.");
				}

				if(dgg.radius==ISEA_SCALE) ret.Append(" +rescale");

				return ret.ToString();
			}
		}

		#region Stuff
		struct hex
		{
			public bool iso;
			public int x, y, z;
		}

		// y *must* be positive down as the xy /iso conversion assumes this
		static void hex_xy(ref hex h)
		{
			if(!h.iso) return;
			if(h.x>=0) h.y=-h.y-(h.x+1)/2;
			else h.y=-h.y-h.x/2; // need to round toward -inf, not toward zero, so x-1

			h.iso=false;

			return;
		}

		static void hex_iso(ref hex h)
		{
			if(h.iso) return;

			if(h.x>=0) h.y=(-h.y-(h.x+1)/2);
			else h.y=(-h.y-(h.x)/2); // need to round toward -inf, not toward zero, so x-1

			h.z=-h.x-h.y;
			h.iso=true;

			return;
		}

		static int hexbin2(double width, double x, double y, out int i, out int j)
		{
			x=x/Math.Cos(30*Math.PI/180.0); // rotated X coord
			y=y-x/2.0; // adjustment for rotated X

			// adjust for actual hexwidth
			x/=width;
			y/=width;

			double z=-x-y;

			double rx=Math.Floor(x+0.5);
			double ry=Math.Floor(y+0.5);
			double rz=Math.Floor(z+0.5);

			int ix=(int)rx;
			int iy=(int)ry;
			int iz=(int)rz;

			int s=ix+iy+iz;

			if(s!=0)
			{
				double abs_dx=Math.Abs(rx-x);
				double abs_dy=Math.Abs(ry-y);
				double abs_dz=Math.Abs(rz-z);

				if(abs_dx>=abs_dy&&abs_dx>=abs_dz) ix-=s;
				else if(abs_dy>=abs_dx&&abs_dy>=abs_dz) iy-=s;
				else iz-=s;
			}

			hex h;
			h.x=ix;
			h.y=iy;
			h.z=iz;
			h.iso=true;

			hex_xy(ref h);
			i=h.x;
			j=h.y;
			return ix*100+iy;
		}

		protected enum isea_poly
		{
			ISEA_NONE,
			ISEA_ICOSAHEDRON=20
		}

		protected enum isea_topology
		{
			ISEA_HEXAGON=6,
			ISEA_TRIANGLE=3,
			ISEA_DIAMOND=4
		}

		protected enum isea_address_form
		{
			ISEA_GEO,
			ISEA_Q2DI,
			ISEA_SEQNUM,
			ISEA_INTERLEAVE,
			ISEA_PLANE,
			ISEA_Q2DD,
			ISEA_PROJTRI,
			ISEA_VERTEX2DD,
			ISEA_HEX
		}

		protected struct isea_dgg
		{
			public isea_poly polyhedron; // ignored, icosahedron
			public double o_lat, o_lon, o_az; // orientation, radians
			//public bool pole; // true if standard snyder
			public isea_topology topology; // ignored, hexagon
			public int aperture; // valid values depend on partitioning method
			public int resolution;
			public double radius; // radius of the earth in meters, ignored 1.0
			public isea_address_form output;
			public int triangle; // triangle of last transformed point
			public int quad; // quad of last transformed point
			public ulong serial;
		}

		struct isea_pt
		{
			public double x, y;
		}

		struct isea_geo
		{
			public double lon, lat;
		}

		//struct isea_address
		//{
		//    public isea_address_form type;
		//    public int number;
		//    public double x, y; // or i, j or lon, lat depending on type
		//}

		enum snyder_polyhedron
		{
			SNYDER_POLY_HEXAGON=0,
			SNYDER_POLY_PENTAGON=1,
			SNYDER_POLY_TETRAHEDRON=2,
			SNYDER_POLY_CUBE=3,
			SNYDER_POLY_OCTAHEDRON=4,
			SNYDER_POLY_DODECAHEDRON=5,
			SNYDER_POLY_ICOSAHEDRON=6
		}

		struct snyder_constants
		{
			public double g, G, theta, ea_w, ea_a, ea_b, g_w, g_a, g_b;
		}

		// TODO put these in radians to avoid a later conversion
		static snyder_constants[] constants=new snyder_constants[]
		{
			new snyder_constants {g=23.80018260, G=62.15458023, theta=60.0, ea_w=3.75, ea_a=1.033, ea_b=0.968, g_w=5.09, g_a=1.195, g_b=1.0},
			new snyder_constants {g=20.07675127, G=55.69063953, theta=54.0, ea_w=2.65, ea_a=1.030, ea_b=0.983, g_w=3.59, g_a=1.141, g_b=1.027},
			new snyder_constants {g=0.0, G=0.0, theta=0.0, ea_w=0.0, ea_a=0.0, ea_b=0.0, g_w=0.0, g_a=0.0, g_b=0.0},
			new snyder_constants {g=0.0, G=0.0, theta=0.0, ea_w=0.0, ea_a=0.0, ea_b=0.0, g_w=0.0, g_a=0.0, g_b=0.0},
			new snyder_constants {g=0.0, G=0.0, theta=0.0, ea_w=0.0, ea_a=0.0, ea_b=0.0, g_w=0.0, g_a=0.0, g_b=0.0},
			new snyder_constants {g=0.0, G=0.0, theta=0.0, ea_w=0.0, ea_a=0.0, ea_b=0.0, g_w=0.0, g_a=0.0, g_b=0.0},
			new snyder_constants {g=37.37736814, G=36.0, theta=30.0, ea_w=17.27, ea_a=1.163, ea_b=0.860, g_w=13.14, g_a=1.584, g_b=1.0},
		};

		const double E=52.62263186;
		const double F=10.81231696;

		const double DEG60=1.04719755119659774614;
		const double DEG120=2.09439510239319549229;
		const double DEG72=1.25663706143591729537;
		const double DEG90=1.57079632679489661922;
		const double DEG144=2.51327412287183459075;
		const double DEG36=0.62831853071795864768;
		const double DEG108=1.88495559215387594306;
		const double DEG180=Math.PI;
		const double ISEA_SCALE=0.8301572857837594396028083; // sqrt(5)/M_PI
		const double V_LAT=0.46364760899944494524; // 26.565051177 degrees
		const double RAD2DEG=(180.0/Math.PI);
		const double DEG2RAD=(Math.PI/180.0);

		static isea_geo[] vertex=new isea_geo[]
		{
			new isea_geo {lon=0.0, lat=DEG90},
			new isea_geo {lon=DEG180, lat=V_LAT},
			new isea_geo {lon=-DEG108, lat=V_LAT},
			new isea_geo {lon=-DEG36, lat=V_LAT},
			new isea_geo {lon=DEG36, lat=V_LAT},
			new isea_geo {lon=DEG108, lat=V_LAT},
			new isea_geo {lon=-DEG144, lat=-V_LAT},
			new isea_geo {lon=-DEG72, lat=-V_LAT},
			new isea_geo {lon=0.0, lat=-V_LAT},
			new isea_geo {lon=DEG72, lat=-V_LAT},
			new isea_geo {lon=DEG144, lat=-V_LAT},
			new isea_geo {lon=0.0, lat=-DEG90}
		};

		// TODO make an isea_pt array of the vertices as well
		static int[] tri_v1= { 0, 0, 0, 0, 0, 0, 6, 7, 8, 9, 10, 2, 3, 4, 5, 1, 11, 11, 11, 11, 11 };

		const double E_RAD=0.91843818702186776133; // 52.62263186
		const double F_RAD=0.18871053072122403508; // 10.81231696

		// triangle Centers
		static isea_geo[] icostriangles=new isea_geo[]
		{
			new isea_geo{lon=0.0, lat=0.0},
			new isea_geo{lon=-DEG144, lat=E_RAD},
			new isea_geo{lon=-DEG72, lat=E_RAD},
			new isea_geo{lon=0.0, lat=E_RAD},
			new isea_geo{lon=DEG72, lat=E_RAD},
			new isea_geo{lon=DEG144, lat=E_RAD},
			new isea_geo{lon=-DEG144, lat=F_RAD},
			new isea_geo{lon=-DEG72, lat=F_RAD},
			new isea_geo{lon=0.0, lat=F_RAD},
			new isea_geo{lon=DEG72, lat=F_RAD},
			new isea_geo{lon=DEG144, lat=F_RAD},
			new isea_geo{lon=-DEG108, lat=-F_RAD},
			new isea_geo{lon=-DEG36, lat=-F_RAD},
			new isea_geo{lon=DEG36, lat=-F_RAD},
			new isea_geo{lon=DEG108, lat=-F_RAD},
			new isea_geo{lon=DEG180, lat=-F_RAD},
			new isea_geo{lon=-DEG108, lat=-E_RAD},
			new isea_geo{lon=-DEG36, lat=-E_RAD},
			new isea_geo{lon=DEG36, lat=-E_RAD},
			new isea_geo{lon=DEG108, lat=-E_RAD},
			new isea_geo{lon=DEG180, lat=-E_RAD},
		};

		static double az_adjustment(int triangle)
		{
			double adj;

			isea_geo v=vertex[tri_v1[triangle]];
			isea_geo c=icostriangles[triangle];

			// TODO looks like the adjustment is always either 0 or 180
			// at least if you pick your vertex carefully
			adj=Math.Atan2(Math.Cos(v.lat)*Math.Sin(v.lon-c.lon),
				Math.Cos(c.lat)*Math.Sin(v.lat)-Math.Sin(c.lat)*Math.Cos(v.lat)*Math.Cos(v.lon-c.lon));

			return adj;
		}

		const double TABLE_G=0.6615845383; // R tan(g) sin(60)
		const double TABLE_H=0.1909830056; // H = 0.25 R tan g =
		const double RPRIME=0.91038328153090290025;

		static isea_pt isea_triangle_xy(int triangle)
		{
			isea_pt c;
			double Rprime=0.91038328153090290025;

			triangle=(triangle-1)%20;

			c.x=TABLE_G*((triangle%5)-2)*2.0;
			if(triangle>9) c.x+=TABLE_G;

			switch(triangle/5)
			{
				case 0: c.y=5.0*TABLE_H; break;
				case 1: c.y=TABLE_H; break;
				case 2: c.y=-TABLE_H; break;
				case 3: c.y=-5.0*TABLE_H; break;
				default:
					// should be impossible
					throw new Exception();
			}

			c.x*=Rprime;
			c.y*=Rprime;

			return c;
		}

		// snyder eq 14
		static double sph_azimuth(double f_lon, double f_lat, double t_lon, double t_lat)
		{
			return Math.Atan2(Math.Cos(t_lat)*Math.Sin(t_lon-f_lon), Math.Cos(f_lat)*Math.Sin(t_lat)-Math.Sin(f_lat)*Math.Cos(t_lat)*Math.Cos(t_lon-f_lon));
		}

		// coord needs to be in radians
		static int isea_snyder_forward(isea_geo ll, out isea_pt @out)
		{
			// TODO by locality of reference, start by trying the same triangle
			// as last time

			// TODO put these constants in as radians to begin with
			snyder_constants c=constants[(int)snyder_polyhedron.SNYDER_POLY_ICOSAHEDRON];

			// plane angle between radius vector to center and adjacent edge of
			// plane polygon
			double theta=c.theta*DEG2RAD;

			// spherical distance from center of polygon face to any of its
			// vertexes on the globe
			double g=c.g*DEG2RAD;

			// spherical angle between radius vector to center and adjacent edge
			// of spherical polygon on the globe
			double G=c.G*DEG2RAD;

			for(int i=1; i<=20; i++)
			{
				double z;
				isea_geo center;

				center=icostriangles[i];

				// step 1
				z=Math.Acos(Math.Sin(center.lat)*Math.Sin(ll.lat)+Math.Cos(center.lat)*Math.Cos(ll.lat)*Math.Cos(ll.lon-center.lon));

				// not on this triangle
				if(z>g+0.000005) continue; // TODO DBL_EPSILON

				double Az=sph_azimuth(center.lon, center.lat, ll.lon, ll.lat);

				// step 2

				// This calculates "some" vertex coordinate
				double az_offset=az_adjustment(i);

				Az-=az_offset;

				// TODO I don't know why we do this. It's not in snyder
				// maybe because we should have picked a better vertex
				if(Az<0.0) Az+=2.0*Math.PI;

				// adjust Az for the point to fall within the range of 0 to
				// 2(90 - theta) or 60 degrees for the hexagon, by
				// and therefore 120 degrees for the triangle
				// of the icosahedron
				// subtracting or adding multiples of 60 degrees to Az and
				// recording the amount of adjustment
				int Az_adjust_multiples=0; // how many multiples of 60 degrees we adjust the azimuth

				while(Az<0.0)
				{
					Az+=DEG120;
					Az_adjust_multiples--;
				}

				while(Az>DEG120+double.Epsilon)
				{
					Az-=DEG120;
					Az_adjust_multiples++;
				}

				// step 3
				double cot_theta=1.0/Math.Tan(theta);
				double tan_g=Math.Tan(g);	// TODO this is a constant

				// Calculate q from eq 9.
				// TODO cot_theta is cot(30)
				double q=Math.Atan2(tan_g, Math.Cos(Az)+Math.Sin(Az)*cot_theta);

				// not in this triangle
				if(z>q+0.000005) continue;
				// step 4

				// Apply equations 5-8 and 10-12 in order

				// eq 5
				// Rprime = 0.9449322893 * R;
				// R' in the paper is for the truncated
				double Rprime=0.91038328153090290025;

				// eq 6
				double H=Math.Acos(Math.Sin(Az)*Math.Sin(G)*Math.Cos(g)-Math.Cos(Az)*Math.Cos(G));

				// eq 7
				// Ag = (Az + G + H - DEG180) * M_PI * R * R / DEG180;
				double Ag=Az+G+H-DEG180;

				// eq 8
				double Azprime=Math.Atan2(2.0*Ag, Rprime*Rprime*tan_g*tan_g-2.0*Ag*cot_theta);

				// eq 10
				// cot(theta) = 1.73205080756887729355
				double dprime=Rprime*tan_g/(Math.Cos(Azprime)+Math.Sin(Azprime)*cot_theta);

				// eq 11
				double f=dprime/(2.0*Rprime*Math.Sin(q/2.0));

				// eq 12
				double rho=2.0*Rprime*f*Math.Sin(z/2.0);

				// add back the same 60 degree multiple adjustment from step
				// 2 to Azprime
				Azprime+=DEG120*Az_adjust_multiples;

				// calculate rectangular coordinates
				double x=rho*Math.Sin(Azprime);
				double y=rho*Math.Cos(Azprime);

				// TODO
				// translate coordinates to the origin for the particular
				// hexagon on the flattened polyhedral map plot
				@out.x=x;
				@out.y=y;

				return i;
			}

			// should be impossible, this implies that the coordinate is not on
			// any triangle

			//fprintf(stderr, "impossible transform: %f %f is not on any triangle\n",
			//    ll.lon*RAD2DEG, ll.lat*RAD2DEG);

			throw new Exception();
		}

		// return the new coordinates of any point in orginal coordinate system.
		// Define a point (newNPold) in orginal coordinate system as the North Pole in
		// new coordinate system, and the great circle connect the original and new
		// North Pole as the lon0 longitude in new coordinate system, given any point
		// in orginal coordinate system, this function return the new coordinates.

		const double PRECISION=0.0000000000005;

		// formula from Snyder, Map Projections: A working manual, p31 */
		//
		// old north pole at np in new coordinates
		// could be simplified a bit with fewer intermediates
		//
		// TODO take a result pointer
		static isea_geo snyder_ctran(isea_geo np, isea_geo pt)
		{
			double phi=pt.lat;
			double lambda=pt.lon;
			double alpha=np.lat;
			double beta=np.lon;
			double lambda0=beta;

			double cos_p=Math.Cos(phi);
			double sin_a=Math.Sin(alpha);

			// mpawm 5-7
			double sin_phip=sin_a*Math.Sin(phi)-Math.Cos(alpha)*cos_p*Math.Cos(lambda-lambda0);

			// mpawm 5-8b

			// lambda prime minus beta
			// use the two argument form so we end up in the right quadrant
			double lp_b=Math.Atan2(cos_p*Math.Sin(lambda-lambda0), (sin_a*cos_p*Math.Cos(lambda-lambda0)+Math.Cos(alpha)*Math.Sin(phi)));

			double lambdap=lp_b+beta;

			// normalize longitude
			// TODO can we just do a modulus ?
			lambdap=Math.IEEERemainder(lambdap, 2*Math.PI);
			while(lambdap>Math.PI) lambdap-=2*Math.PI;
			while(lambdap<-Math.PI) lambdap+=2*Math.PI;

			double phip=Math.Asin(sin_phip);

			isea_geo npt;

			npt.lat=phip;
			npt.lon=lambdap;

			return npt;
		}

		static isea_geo isea_ctran(ref isea_geo np, isea_geo pt, double lon0)
		{
			isea_geo npt;

			np.lon+=Math.PI;
			npt=snyder_ctran(np, pt);
			np.lon-=Math.PI;

			npt.lon-=(Math.PI-lon0+np.lon);

			// snyder is down tri 3, isea is along side of tri1 from vertex 0 to
			// vertex 1 these are 180 degrees apart
			npt.lon+=Math.PI;

			// normalize longitude
			npt.lon=Math.IEEERemainder(npt.lon, 2*Math.PI);
			while(npt.lon>Math.PI) npt.lon-=2*Math.PI;
			while(npt.lon<-Math.PI) npt.lon+=2*Math.PI;

			return npt;
		}

		// in radians
		const double ISEA_STD_LAT=1.01722196792335072101;
		const double ISEA_STD_LON=0.19634954084936207740;

		// fuller's at 5.2454 west, 2.3009 N, adjacent at 7.46658 deg

		static void isea_grid_init(ref isea_dgg g)
		{
			g.polyhedron=isea_poly.ISEA_ICOSAHEDRON;
			g.o_lat=ISEA_STD_LAT;
			g.o_lon=ISEA_STD_LON;
			g.o_az=0.0;
			g.aperture=4;
			g.resolution=6;
			g.radius=1.0;
			g.topology=isea_topology.ISEA_HEXAGON;
		}

		static void isea_orient_isea(ref isea_dgg g)
		{
			g.o_lat=ISEA_STD_LAT;
			g.o_lon=ISEA_STD_LON;
			g.o_az=0.0;
		}

		static void isea_orient_pole(ref isea_dgg g)
		{
			g.o_lat=Math.PI/2.0;
			g.o_lon=0.0;
			g.o_az=0;
		}

		static int isea_transform(ref isea_dgg g, isea_geo @in, out isea_pt @out)
		{
			isea_geo pole;
			pole.lat=g.o_lat;
			pole.lon=g.o_lon;

			isea_geo i=isea_ctran(ref pole, @in, g.o_az);

			int tri=isea_snyder_forward(i, out @out);
			@out.x*=g.radius;
			@out.y*=g.radius;
			g.triangle=tri;

			return tri;
		}

		static void isea_rotate(ref isea_pt pt, double degrees)
		{
			double rad=-degrees*Math.PI/180.0;
			while(rad>=2.0*Math.PI) rad-=2.0*Math.PI;
			while(rad<=-2.0*Math.PI) rad+=2.0*Math.PI;

			double x=pt.x*Math.Cos(rad)+pt.y*Math.Sin(rad);
			double y=-pt.x*Math.Sin(rad)+pt.y*Math.Cos(rad);

			pt.x=x;
			pt.y=y;
		}

		static int isea_tri_plane(int tri, ref isea_pt pt, double radius)
		{
			isea_pt tc; // center of triangle

			if(((tri-1)/5)%2==1) isea_rotate(ref pt, 180.0);

			tc=isea_triangle_xy(tri);
			tc.x*=radius;
			tc.y*=radius;
			pt.x+=tc.x;
			pt.y+=tc.y;

			return tri;
		}

		// convert projected triangle coords to quad xy coords, return quad number
		static int isea_ptdd(int tri, ref isea_pt pt)
		{
			bool downtri=(((tri-1)/5)%2==1);
			int quad=((tri-1)%5)+((tri-1)/10)*5+1;

			isea_rotate(ref pt, downtri?240.0:60.0);
			if(downtri)
			{
				pt.x+=0.5;
				// pt->y += cos(30.0 * Math.PI / 180.0);
				pt.y+=0.86602540378443864672;
			}
			return quad;
		}

		static int isea_dddi_ap3odd(ref isea_dgg g, int quad, isea_pt pt, out isea_pt di)
		{
			// This is the number of hexes from apex to base of a triangle
			double sidelength=(Math.Pow(2.0, g.resolution)+1.0)/2.0; // in hexes

			// apex to base is cos(30deg)
			double hexwidth=Math.Cos(Math.PI/6.0)/sidelength;

			// TODO I think sidelength is always x.5, so
			// (int)sidelength * 2 + 1 might be just as good
			int maxcoord=(int)(sidelength*2.0+0.5);

			isea_pt v=pt;
			hex h;
			h.z=0;
			hexbin2(hexwidth, v.x, v.y, out h.x, out h.y);
			h.iso=false;
			hex_iso(ref h);

			int d=h.x-h.z;
			int i=h.x+h.y+h.y;

			// you want to test for max coords for the next quad in the same
			// "row" first to get the case where both are max
			if(quad<=5)
			{
				if(d==0&&i==maxcoord)
				{
					// north pole
					quad=0;
					d=0;
					i=0;
				}
				else if(i==maxcoord)
				{
					// upper right in next quad
					quad+=1;
					if(quad==6) quad=1;
					i=maxcoord-d;
					d=0;
				}
				else if(d==maxcoord)
				{
					// lower right in quad to lower right
					quad+=5;
					d=0;
				}
			}
			else if(quad>=6)
			{
				if(i==0&&d==maxcoord)
				{
					// south pole
					quad=11;
					d=0;
					i=0;
				}
				else if(d==maxcoord)
				{
					// lower right in next quad
					quad+=1;
					if(quad==11) quad=6;
					d=maxcoord-i;
					i=0;
				}
				else if(i==maxcoord)
				{
					// upper right in quad to upper right
					quad=(quad-4)%5;
					i=0;
				}
			}

			di.x=d;
			di.y=i;

			g.quad=quad;
			return quad;
		}

		static int isea_dddi(ref isea_dgg g, int quad, isea_pt pt, out isea_pt di)
		{
			if(g.aperture==3&&g.resolution%2!=0)
				return isea_dddi_ap3odd(ref g, quad, pt, out di);

			// todo might want to do this as an iterated loop
			int sidelength;	// in hexes
			if(g.aperture>0) sidelength=(int)(Math.Pow(g.aperture, g.resolution/2.0)+0.5);
			else sidelength=g.resolution;

			double hexwidth=1.0/sidelength;

			isea_pt v=pt;
			hex h;
			h.z=0;
			isea_rotate(ref v, -30.0);
			hexbin2(hexwidth, v.x, v.y, out h.x, out h.y);
			h.iso=false;
			hex_iso(ref h);

			// we may actually be on another quad
			if(quad<=5)
			{
				if(h.x==0&&h.z==-sidelength)
				{
					// north pole
					quad=0;
					h.z=0;
					h.y=0;
					h.x=0;
				}
				else if(h.z==-sidelength)
				{
					quad=quad+1;
					if(quad==6) quad=1;
					h.y=sidelength-h.x;
					h.z=h.x-sidelength;
					h.x=0;
				}
				else if(h.x==sidelength)
				{
					quad+=5;
					h.y=-h.z;
					h.x=0;
				}
			}
			else if(quad>=6)
			{
				if(h.z==0&&h.x==sidelength)
				{
					// south pole
					quad=11;
					h.x=0;
					h.y=0;
					h.z=0;
				}
				else if(h.x==sidelength)
				{
					quad=quad+1;
					if(quad==11) quad=6;
					h.x=h.y+sidelength;
					h.y=0;
					h.z=-h.x;
				}
				else if(h.y==-sidelength)
				{
					quad-=4;
					h.y=0;
					h.z=-h.x;
				}
			}
			di.x=h.x;
			di.y=-h.z;

			g.quad=quad;
			return quad;
		}

		static int isea_ptdi(ref isea_dgg g, int tri, isea_pt pt, out isea_pt di)
		{
			isea_pt v=pt;
			int quad=isea_ptdd(tri, ref v);
			quad=isea_dddi(ref g, quad, v, out di);
			return quad;
		}

		// q2di to seqnum
		static int isea_disn(ref isea_dgg g, int quad, isea_pt di)
		{
			if(quad==0)
			{
				g.serial=1;
				return 1;
			}

			// hexes in a quad
			int hexes=(int)(Math.Pow(g.aperture, g.resolution)+0.5);
			if(quad==11)
			{
				g.serial=(ulong)(1+10*hexes+1);
				return (int)g.serial;
			}

			int sn;
			if(g.aperture==3&&g.resolution%2==1)
			{
				int height=(int)(Math.Pow(g.aperture, (g.resolution-1)/2.0));
				sn=((int)di.x)*height;
				sn+=((int)di.y)/height;
				sn+=(quad-1)*hexes;
				sn+=2;
			}
			else
			{
				int sidelength=(int)(Math.Pow(g.aperture, g.resolution/2.0)+0.5);
				sn=(int)((quad-1)*hexes+sidelength*di.x+di.y+2);
			}

			g.serial=(ulong)sn;
			return sn;
		}

		// TODO just encode the quad in the d or i coordinate
		// quad is 0-11, which can be four bits.
		// d' = d << 4 + q, d = d' >> 4, q = d' & 0xf

		// convert a q2di to global hex coord
		static int isea_hex(ref isea_dgg g, int tri, isea_pt pt, out isea_pt hex)
		{
			isea_pt v;
			int quad=isea_ptdi(ref g, tri, pt, out v);

			hex.x=((int)v.x<<4)+quad;
			hex.y=v.y;

			return 1;

// silence the compiler
#if false
			double d=v.x;
			double i=v.y;

			// Aperture 3 odd resolutions
			if(g.aperture==3&&g.resolution%2!=0)
			{
				int offset=(int)(Math.Pow(3.0, g.resolution-1)+0.5);

				d+=offset*((g.quad-1)%5);
				i+=offset*((g.quad-1)%5);

				if(quad==0)
				{
					d=0;
					i=offset;
				}
				else if(quad==11)
				{
					d=2*offset;
					i=0;
				}
				else if(quad>5)
				{
					d+=offset;
				}

				double x=(2*d-i)/3;
				double y=(2*i-d)/3;

				hex.x=x+offset/3;
				hex.y=y+2*offset/3;

				return 1;
			}

			// aperture 3 even resolutions and aperture 4
			double sidelength=(int)(Math.Pow(g.aperture, g.resolution/2.0)+0.5);
			if(g.quad==0)
			{
				hex.x=0;
				hex.y=sidelength;
			}
			else if(g.quad==11)
			{
				hex.x=sidelength*2;
				hex.y=0;
			}
			else
			{
				hex.x=d+sidelength*((g.quad-1)%5);
				if(g.quad>5) hex.x+=sidelength;
				hex.y=i+sidelength*((g.quad-1)%5);
			}

			return 1;
#endif
		}

		static isea_pt isea_forward(isea_dgg g, isea_geo @in)
		{
			isea_pt @out, coord;

			int tri=isea_transform(ref g, @in, out @out);

			if(g.output==isea_address_form.ISEA_PLANE)
			{
				isea_tri_plane(tri, ref @out, g.radius);
				return @out;
			}

			// convert to isea standard triangle size
			@out.x=@out.x/g.radius*ISEA_SCALE;
			@out.y=@out.y/g.radius*ISEA_SCALE;
			@out.x+=0.5;
			@out.y+=2.0*.14433756729740644112;

			switch(g.output)
			{
				case isea_address_form.ISEA_PROJTRI:
					// nothing to do, already in projected triangle
					break;
				case isea_address_form.ISEA_VERTEX2DD:
					g.quad=isea_ptdd(tri, ref @out);
					break;
				case isea_address_form.ISEA_Q2DD:
					// Same as above, we just don't print as much
					g.quad=isea_ptdd(tri, ref @out);
					break;
				case isea_address_form.ISEA_Q2DI:
					g.quad=isea_ptdi(ref g, tri, @out, out coord);
					return coord;
				case isea_address_form.ISEA_SEQNUM:
					isea_ptdi(ref g, tri, @out, out coord);
					// disn will set g->serial
					isea_disn(ref g, g.quad, coord);
					return coord;
				case isea_address_form.ISEA_HEX:
					isea_hex(ref g, tri, @out, out coord);
					return coord;
			}

			return @out;
		}
		#endregion

		// Proj 4 integration code follows

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			isea_geo @in;
			@in.lon=lp.lam;
			@in.lat=lp.phi;

			isea_pt @out=isea_forward(dgg, @in);

			xy.x=@out.x;
			xy.y=@out.y;

			return xy;
		}

		public override PJ Init()
		{
			fwd=s_forward;
			isea_grid_init(ref dgg);

			dgg.output=isea_address_form.ISEA_PLANE;
			//dgg.radius = a; // otherwise defaults to 1
			// calling library will scale, I think

			string opt=Proj.pj_param_s(ctx, parameters, "orient");

			if(opt!=null&&opt!="")
			{
				if(opt=="isea") isea_orient_isea(ref dgg);
				else if(opt=="pole") isea_orient_pole(ref dgg);
				else
				{
					Proj.pj_ctx_set_errno(ctx, -34);
					return null;
				}
			}

			if(Proj.pj_param_t(ctx, parameters, "azi"))
				dgg.o_az=Proj.pj_param_r(ctx, parameters, "azi");

			if(Proj.pj_param_t(ctx, parameters, "lon_0"))
				dgg.o_lon=Proj.pj_param_r(ctx, parameters, "lon_0");

			if(Proj.pj_param_t(ctx, parameters, "lat_0"))
				dgg.o_lat=Proj.pj_param_r(ctx, parameters, "lat_0");

			if(Proj.pj_param_t(ctx, parameters, "aperture"))
				dgg.aperture=Proj.pj_param_i(ctx, parameters, "aperture");
			else dgg.aperture=3;

			if(Proj.pj_param_t(ctx, parameters, "resolution"))
				dgg.resolution=Proj.pj_param_i(ctx, parameters, "resolution");
			else dgg.resolution=4;

			opt=Proj.pj_param_s(ctx, parameters, "mode");
			if(opt!=null&&opt!="")
			{
				if(opt=="plane") dgg.output=isea_address_form.ISEA_PLANE;
				else if(opt=="di") dgg.output=isea_address_form.ISEA_Q2DI;
				else if(opt=="dd") dgg.output=isea_address_form.ISEA_Q2DD;
				else if(opt=="hex") dgg.output=isea_address_form.ISEA_HEX;
				else
				{
					// TODO verify error code. Possibly eliminate magic
					Proj.pj_ctx_set_errno(ctx, -34);
					return null;
				}
			}

			if(Proj.pj_param_t(ctx, parameters, "rescale"))
				dgg.radius=ISEA_SCALE;

			return this;
		}
	}
}

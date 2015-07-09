namespace Free.Ports.Proj4.Gridshift
{
	public class CTABLE
	{
		public string id;	// ascii info
		public LP ll;		// lower left corner coordinates
		public LP del;		// size of cells
		public ILP lim;		// limits of conversion matrix
		public LP[] cvs;	// conversion matrix

		public CTABLE Clone()
		{
			CTABLE ret=new CTABLE();
			ret.id=id;
			ret.ll=ll;
			ret.del=del;
			ret.lim=lim;
			ret.cvs=cvs;
			return ret;
		}
	}
}

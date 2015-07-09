namespace Free.Ports.Proj4.FactorsDeriv
{
	public delegate void void_LP_PJ_FACTORS(LP lp, ref FACTORS fac);

	public interface IFactors
	{
		void_LP_PJ_FACTORS spc { get; }
	}
}

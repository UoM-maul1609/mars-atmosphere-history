"""Run: python -m unittest discover -s tests -v (from repo root)."""
import os
import sys
import unittest
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import MarsSputtering_Calc6_7 as sput
import isotopic_data as iso

class PaperEquationTests(unittest.TestCase):
    def test_table2_nodes(self):
        for euv,z_km,t_k in [(1,80,200),(3,125,300),(6,175,500),(100,600,1200)]:
            z,t=sput.upper_atmosphere(euv)
            self.assertAlmostEqual(z/1000,z_km)
            self.assertAlmostEqual(t,t_k)
    def test_reference_masses(self):
        self.assertEqual(sput.amus['C'],44)
        self.assertEqual(sput.amus['N'],28)
        self.assertAlmostEqual(sput.R_diffij(0.,sput.g_mars,80e3,200),1)
    def test_equal_factors_keep_isotopes_unchanged(self):
        ratios=[0.2,0.3,0.4]
        out=iso._escape_one(ratios,100,15,[1,1,1])
        np.testing.assert_allclose(out,ratios,rtol=1e-12)
    def test_mixing_preserves_isotope_populations(self):
        old=[0.2,0.4]
        source=[0.5,0.1]
        out=iso._mix_one(old,source,75.,25.)
        expected=75*iso._weights(old)+25*iso._weights(source)
        np.testing.assert_allclose(iso._weights(out),expected/100,rtol=1e-12)
    def test_heavy_enrichment_with_escape(self):
        out=iso._escape_one([0.003],100,10,[0.7])
        self.assertGreater(out[0],0.003)
    def test_n_photochemical_factors(self):
        r=0.003
        sput_only=iso.escape_fr([[r]], [100.], [10.], [0.], [[0.9]])[0][0]
        photo_only=iso.escape_fr([[r]], [100.], [0.], [10.], [[0.9]])[0][0]
        self.assertGreater(photo_only,sput_only)
    def test_diffusion_heavy_isotope_escape_suppressed(self):
        # Heavier isotopes have a smaller exobase/homopause ratio.
        ref = [14,20,36,84,130]
        heavy = [[15],[22],[38],[86],[136]]
        for euv in (1.,3.,6.,100.):
            factors = sput.calculate_rdiffs(ref,heavy,euv)
            for f in factors:
                self.assertGreater(f[0],0.)
                self.assertLess(f[0],1.)
        self.assertAlmostEqual(sput.calculate_rdiffs([14],[[14]])[0][0],1.)

    def test_sputtering_yield_co2_reference(self):
        # Pure CO2: Rdiff=alpha=1, so carbon/CO2 loss matches FCO2.
        res=sput.F_i_sp(1e25,np.array([1e20]),1e20,['C'],euv=1.)
        self.assertAlmostEqual(res[0]/1e25,1.,places=12)

    def test_reference_basis_mixing_weights_reference_isotope(self):
        # Reference-isotope ledgers: 130Xe-like counts add directly.
        atm=[16.4, 2.3]; src=[6.4, 2.0]
        out=iso._mix_one(atm,src,30.,10.,reference_basis=True)
        expected=(30*np.array(atm)+10*np.array(src))/40.
        np.testing.assert_allclose(out,expected,rtol=1e-12)

    def test_reference_basis_escape_removes_loss_of_reference(self):
        # Equal factors: ratios unchanged; reference loses exactly `loss`.
        out=iso._escape_one([0.5,2.0],100.,10.,[1.,1.],reference_basis=True)
        np.testing.assert_allclose(out,[0.5,2.0],rtol=1e-12)

    def test_idp_rates_are_reference_isotope_moles(self):
        import IDPs
        rate=IDPs.get_mole(['Xe'])[0]
        f_ref,mean_mass=IDPs._reference_fraction_and_mean_mass('Xe')
        self.assertAlmostEqual(rate, 7.33e9/mean_mass*f_ref)
        self.assertLess(f_ref,0.05)  # 130Xe is ~4.4% of solar-like Xe

    def test_co2_yield_matches_table1(self):
        self.assertEqual(sput.Y_CO2, 0.7)
        self.assertEqual(sput.yields['C'], 0.7)

    # ---- discrete impacts -------------------------------------------------
    def test_truncated_sizes_within_bounds_and_mean(self):
        import impactor_model as im
        # caps below and above the 300 km slope break
        for dmax in (100e3, 500e3):
            d,m=im.sample_sizes_truncated(np.random.RandomState(0).rand(2000000),dmax=dmax)
            self.assertGreaterEqual(d.min(),im.dmin); self.assertLessEqual(d.max(),dmax)
            self.assertAlmostEqual(m.mean()/im.mean_impactor_mass(dmax=dmax),1.,delta=0.03)

    def test_impact_count_matches_kurokawa_fig3a(self):
        # Fig. 3a: ~2.3e5 impacts with D >= 3 km for 2e21 kg.
        import impactor_model as im
        n=im.number_of_impacts(dmax=im.budget_diameter())
        self.assertTrue(1.8e5 < n < 3.2e5, n)

    def test_impact_times_follow_weights(self):
        import impactor_model as im
        st,d,m=im.sample_impact_history(np.random.RandomState(1),[1.,0.,3.])
        self.assertTrue(np.all(np.diff(st)>=0))
        self.assertNotIn(1,st)
        frac=np.mean(st==2); self.assertAlmostEqual(frac,0.75,delta=0.01)

    def test_run_is_reproducible_for_a_seed(self):
        import main
        a=main.run_model(3,-1)[1][-1]; b=main.run_model(3,-1)[1][-1]
        np.testing.assert_array_equal(a,b)

    # ---- nitrogen and open-system escape ---------------------------------
    def test_open_system_steady_state_and_conservation(self):
        # Fast loss with supply relaxes to s/k with ratio = source / R.
        r,inv,rem=iso.evolve_open_system([0.2],100.,[1.0],[[0.8]],[(50.,[0.2])],1e3,True)
        self.assertAlmostEqual(r[0],0.25,places=10); self.assertAlmostEqual(inv,0.05,places=10)
        # Mass balance: initial + supplied = final + removed (total basis).
        r,inv,rem=iso.evolve_open_system([0.01],100.,[0.3,0.2],[[0.9],[0.6]],[(40.,[0.004])],2.,False)
        self.assertAlmostEqual(100.+40.,inv+sum(rem),places=9)
        # No loss reduces to exact mixing.
        np.testing.assert_allclose(iso.evolve_open_system([0.3],70.,[0.],[[1.]],[(30.,[0.1])],1.,True)[0],
                                   iso._mix_one([0.3],[0.1],70.,30.,True))

    def test_photochemical_N_uses_bulk_ratio(self):
        # At the present bulk N2/CO2 (0.0264) the flux equals the tabulated
        # present-day Fox & Dalgarno sum at 1x EUV (phi = 1).
        import photo_chemical as pc
        A=1.0
        f=sum(pc.nitrogen_escape(1.0,0.026390532,1.,1.,A,0.))
        # (digitised phi curves read 1.0026 at the normalisation point)
        for act,total in (('low',2.26e5),('average',5.494e5)):
            pc.FD83_SOLAR_ACTIVITY=act
            f=sum(pc.nitrogen_escape(1.0,0.026390532,1.,1.,A,0.))
            self.assertAlmostEqual(f/1e4/total,1.,delta=0.006)
        pc.FD83_SOLAR_ACTIVITY='low'

    def test_euv_and_sputtering_nodes(self):
        import sputtering_co2 as sc
        self.assertAlmostEqual(float(sc.euv_relative(4.5)),1.,places=6)
        rates=sc.CO2_SPUTTERING_RATES[sc.CO2_SPUTTERING_SOURCE]
        for lvl,f in ((6.,rates[0]),(3.,rates[1])):
            e,F=sc.sputtering_co2_rate(sc.time_of_euv(lvl))
            self.assertAlmostEqual(float(e),lvl,places=6); self.assertAlmostEqual(float(F)/f,1.,places=6)

    def test_lifetimes_match_kurokawa_fig14(self):
        # Kurokawa Fig. 14 at 6 mbar, N2/CO2 = 1e-2 (Ar/CO2 = 1e-6), read off
        # the published contours: N photochemical ~2 Gyr, N sputtering ~2 Gyr,
        # Ar sputtering ~6 Gyr at present; N sputtering ~5 Myr at 3.5 Ga.
        import photo_chemical as pc, sputtering_co2 as sc
        A=4*np.pi*3389.5e3**2; GYR=3.156e16; NA=6.02e23
        def inv(r):
            n=6e-3*1e5*A/3.721/(0.044+r*0.028)*NA; return n, r*n
        e,F=sc.sputtering_co2_rate(4.5)
        co2,n2=inv(1e-2)
        tph=2*n2/sum(pc.nitrogen_escape(co2,n2,e,1.,A,0))/GYR
        tsp=n2/sput.F_i_sp(F,np.array([co2,n2]),co2,['C','N'],euv=e)[1]/GYR
        co2b,ar=inv(1e-6)
        tar=ar/sput.F_i_sp(F,np.array([co2b,ar]),co2b,['C','Ar'],euv=e)[1]/GYR
        for got,want in ((tph,2.),(tsp,2.),(tar,6.)):
            self.assertTrue(0.6 < got/want < 1.6, (got,want))
        e,F=sc.sputtering_co2_rate(4.5-3.5)
        tsp35=n2/sput.F_i_sp(F,np.array([co2,n2]),co2,['C','N'],euv=e)[1]/GYR
        self.assertTrue(0.6 < tsp35/0.005 < 1.6, tsp35)

    def test_fig2_normalisation_and_comet_volatiles(self):
        # Fig. 2 asteroids should come out at ~1 wt% CO2 (the X_gas = 1%
        # case); comets at ~3x that, independent of X_gas.
        Kkg=iso.FIG2_MOL_PER_KG
        ast=iso.abundances1['Asteroids'][1][0]*iso.solar_abundances[0]*Kkg*44.01e-3
        self.assertTrue(0.007 < ast < 0.012, ast)
        self.assertTrue(0.025 < iso.comet_co2_mass_fraction() < 0.033)

    def test_no_default_upper_pressure_cap(self):
        # Kurokawa Fig. 5a: X_gas = 10% gives tens of bar; no 2.75 bar cap.
        import main
        y=main.run_model(0,-1,X_gas=0.1,C_vol=1.,C_Ne_IDP=1.)[1]
        self.assertGreater(y[-1,0], 3.)
        out=main.run_model(0,-1,X_gas=0.1,C_vol=1.,C_Ne_IDP=1.,P_upper_cap=2.75)
        y,A=out[1],out[4]
        pco2=y[:,8]*3.721/A/1e5          # cap applies to CO2, not total
        self.assertLessEqual(pco2[1:].max(), 2.75+1e-9)

    def test_noble_gas_tracers_sum_to_inventory(self):
        import main
        tr={}
        me=main.run_model(14,-1,tracers=tr)[2]
        np.testing.assert_allclose(sum(tr.values()), me[2:6], rtol=1e-8)

    def test_regolith_reservoir_and_xe_ion_escape(self):
        import main
        base=main.run_model(3,-1)
        # zero K and no ion escape reproduce the default exactly
        same=main.run_model(3,-1,K_reg=[0,0,0,0],xe_ion_k=0.)
        np.testing.assert_allclose(same[2], base[2], rtol=1e-12)
        # ion escape fractionates Xe only (from a chondritic start)
        # (checked at 4.1 Ga: with K = 0 little early Xe survives to today)
        esc=main.run_model(3,-1,xe_ion_k=5.)
        ref=main.run_model(3,-1,xe_initial='Volcanic degassing')
        j=np.argmin(abs(esc[0]+4.1e9))
        self.assertGreater(esc[1][j,5], ref[1][j,5] + 0.05)
        self.assertAlmostEqual(esc[1][-1,4], base[1][-1,4], places=10)
        # a regolith reservoir shields early (pre-4.1 Ga) Xe from erosion
        tr0,tr1={},{}
        main.run_model(3,-1,tracers=tr0)
        main.run_model(3,-1,tracers=tr1,K_reg=[0,0,0,30])
        f=lambda tr: tr['pre'][3]/sum(v[3] for v in tr.values())
        self.assertGreater(f(tr1), 2*f(tr0))

if __name__ == '__main__': unittest.main()

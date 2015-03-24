Models with full grid data:
* Gravel, 60ky split time, just simulate East Asians: null_model_gravel_asn_scale_60k.sh
* Simple two-pop model with 60ky split, Ne=10k, no migration: null_model_two_pop_10k_10k_base.sh
* Simple two-pop model with 60ky split, Ne_YRI=10k, Ne_nonAfrican=5k, no migration: null_model_two_pop_10k_5k_base.sh
* Simple two-pop two-wave OOA model with 60ky first wave, 50ky second wave, and joining again at 4kya. Ne_YRI=10k, Ne_nonAfrican=5k, no migration: null_model_two_pop_two_wave_10k_5k_base.sh

Test models:
* EUR, Gravel, larger Ne: test_run_null_simulations.seg.sh
* EUR, base Gravel: test_run_null_simulations.seg.gravel.sh
* ASN, base Gravel: test_run_null_simulations.seg.gravel.asn.sh
* ASN, Gravel, EUR/ASN split at 45kya: test_run_null_simulations.seg.gravel.asn.scale.sh
* EUR: Gravel, EUR/ASN split at 45kya: test_run_null_simulations.seg.gravel.eur.scale.sh

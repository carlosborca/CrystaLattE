import crystalatte
from crystalatte import plugins


def func():
    crystalatte.main(
        cif_input="./Tests/Ammonia/Ammonia.cif",
        cif_output="./Tests/Ammonia/Ammonia.xyz",
        cif_a=3,
        cif_b=3,
        cif_c=3,
        bfs_thresh=1.2,
        uniq_filter="ChSEV",
        nmers_up_to=3,
        r_cut_com=6.5,
        r_cut_monomer=3.5,
        r_cut_dimer=3.6,
        r_cut_trimer=5.7,
        r_cut_tetramer=3.7,
        r_cut_pentamer=6.1,
        cle_run_type=["custom"],
        method="my_method",
        bsse_type=None,
        job_memory=None,
        verbose=2,
        custom_function=plugins.force_fields.example_energy_function
    )
    return


def main():
    func()
    return


if __name__ == "__main__":
    main()

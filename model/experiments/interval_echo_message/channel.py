"""Complete terminal channel on the intact federation."""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate as echo
from model.experiments.interval_echo_space.instrument import local_tables

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SPEC = json.loads((HERE / "protocol.json").read_text())
TOL = SPEC["numeric_tolerance"]

# Whole-file changes reviewed at OPH 0ebb9ba9d6280fec258a7d65e2966c13352dfa8d.
REPLAY_UPDATES = {
    "model/candidate_laws/interval_echo/SPEC.md": "46781970157edaec4b30e55bb673c90102ffcf937aeb6974927e58b08e74dff4",
    "model/candidate_laws/interval_echo/candidate.py": "a318f231c29753da6d0bfaf54ceb83f4410ede7b054b1edc1360620692876944",
    "model/candidate_laws/interval_echo/EVALUATION.json": "4a31cbd583995fa9af119b4e4a4f5545c3656f5f26dfd25de49d6a9440509331",
}

# Additional metadata changes reviewed at main 70e58e24f3be9295926164d38626d58acb87b719.
REPLAY_70E58E24 = {
    "model/candidate_laws/interval_echo/SPEC.md": "c9954a39eaef1de8788e5c4a0e1de1cbd81e8d4dd162566c6aac59467d1af506",
    "model/candidate_laws/interval_echo/EVALUATION.json": "7d95fd5f067232b9f9f3aba664148865ee6773b90ffa52e92f882b06bdae80d9",
}

# The regional conservation, coherent-sweep and coupling evidence changes metadata only.
REPLAY_REGIONAL_CONSENSUS = {
    "model/candidate_laws/interval_echo/candidate.py": "c55b45990fedd5b79ebf5d8e206334b7e90c519e72811b5687149c7f58d9a4a4",
    "model/candidate_laws/interval_echo/EVALUATION.json": "768e20a009f611483bd51a027cc47593a13084e276fb35a561a9cde6953d94cc",
}

# Candidate and evaluation metadata now also record the separate comparison standard.
REPLAY_MODULAR_STANDARD = {
    "model/candidate_laws/interval_echo/candidate.py": "6afe9f5a24481beb0331a7f6930fce6d22c72cb6147d4491aa68b49edeb6107d",
    "model/candidate_laws/interval_echo/SPEC.md": "846dc2cbe46a714dd87d735d3a9d4d21e11bd022da8a6d106b3b351ce0c5c00a",
    "model/candidate_laws/interval_echo/EVALUATION.json": "661a1fbd55b9766acad1b657da24578b3a2ab35cad2016aad129a5d6d165cce6"
}


# Reviewed premise and evidence metadata for the regional gravitational and dark-sector instruments.
REPLAY_GRAVITY_CAP = {
    "model/candidate_laws/interval_echo/candidate.py": "c69c10a59f40c5e71f20854d0fec80faa5a3be16cf4f9232adae884dadb0b59e",
    "model/candidate_laws/interval_echo/EVALUATION.json": "c2817432c2218a97220d56d107f8c81748178f5c5eb8788de665ec291e46920f"
}


# Nested cap, boundary-volume and clock comparisons update the premise/evidence metadata.
REPLAY_NESTED_CAP = {
    "model/candidate_laws/interval_echo/candidate.py": "c666a31015ea2203ac919d55d98ddd5e47fe78076b3d87dd94bed2dd9990a503",
    "model/candidate_laws/interval_echo/EVALUATION.json": "6ad745ce68fa19b1c2fbdc1ea31948403e1d33bc13b8b761440141a031b2661a",
}


# G6 now records the joint clock and load readout on the same history.
REPLAY_JOINT_READ = {
    "model/candidate_laws/interval_echo/EVALUATION.json": "801e9bb3b1ae0ed081c74a5685a6b01c9c41410239206211ade58e1fbccd3d17",
}


# Directional modular charge, material radar and source discrimination extend evaluation metadata.
REPLAY_DIRECTIONAL_CHARGE = {
    "model/candidate_laws/interval_echo/EVALUATION.json": "8e59ded545e23fd3419c224933ca326ee0bc9926a023db87e85457e438ba32f4",
}


# Each reviewed metadata update is matched against the complete law and instrument manifest.
REPLAY_EVALUATIONS = (
    {"model/candidate_laws/interval_echo/EVALUATION.json": "7a1b0f8fef54157548ac1c40c4e107934b766d5f5994b94791506fd77a245530"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "41c5457f88cc152a719a8678a783838c00d6255b62fd9cfbed3379aedb61683c"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "8e50a50244e99f851faf42720b1e72555fc52f2508fb3eef46c23412736266e0"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "3aee1acfeeea03479956235d0050073b21f255b9389bf3dce2fc32384ffd017b"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "5fded60f359db42455f4977d7be03134ea74afa9109dc2bf242de95573a6d128"},
    # Row 4 filing status and station geometry receipts are evaluation metadata only.
    {"model/candidate_laws/interval_echo/EVALUATION.json": "5645b804c85ff2107a73b970b15973a18c81e2c45e2881a16639ac3caa1668af"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "284edfcb6c696e668563c09d7134b6d644c7784ec430955c9b0cecf331589388"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "11cbdf312c1265dd8c46330d55703c2f9e54bf6bfae9e885f8f4bb500e22029a"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "bed52a6bf26421a12a2e424c36270000febf3bfa09b10e9f7490a597b4ebabc0"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "f5e5a7c4dd57c9ec39ebd26f56299708bb6788f8a2fbcfeb6ee4dcddf495f12e"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "c3d796e939d8223d0f59494cbc9eb8344295150a142bfc68f5813151b072b3bd"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "03513c4549fd56c40c2867da9c51ae9476b2e68a2c6657dc497883ab175839ae"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "2414fb6856f0472fdba6de899155f1f66073ac1629617430cf125238e5f00a32"},
    {"model/candidate_laws/interval_echo/EVALUATION.json": "2364858507fdf8615607557a4d4dc191f170b79b2f054155f882170fe7bf384f"},
    # The native message reader replays the complete metadata profile at 5f3df4f03.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "9b08035ebb2796af7227cfda6915b2f354fb1d2d97bdecaa47ecfeada1941eec",
        "model/candidate_laws/interval_echo/candidate.py": "7cd3ab08a927cf89c1323922091f92b505151a693502d46d3d86c4dace2b1b1c",
        "model/candidate_laws/interval_echo/EVALUATION.json": "5f99c88503cf532871010f207af0e77d3e7f48f57d17eb4942bf6be84208acc1",
    },
    # Source-bound orbital, colour and gravity readouts extend candidate and scoring metadata.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "9b08035ebb2796af7227cfda6915b2f354fb1d2d97bdecaa47ecfeada1941eec",
        "model/candidate_laws/interval_echo/candidate.py": "babb7b616a19021a67f7f0d7ab146d7e41c0f2a6923f5b7e9a8578e28916ddef",
        "model/candidate_laws/interval_echo/EVALUATION.json": "ede7beb306d9f46eba035cf94f0b8fc09fe604a92b8d1bf55076b40242a23940",
    },
    # The loop-transport paragraph extends the candidate description only.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "4a6c04fb18a52dca23c33e74dc50d9c2e3d24939e710e5ba1583c8797d7460ff",
        "model/candidate_laws/interval_echo/candidate.py": "babb7b616a19021a67f7f0d7ab146d7e41c0f2a6923f5b7e9a8578e28916ddef",
        "model/candidate_laws/interval_echo/EVALUATION.json": "ede7beb306d9f46eba035cf94f0b8fc09fe604a92b8d1bf55076b40242a23940",
    },
    # Particle identifications and finite measurement reports extend the evaluation metadata.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "79fd7900f8b7dd8a8ba6aac02dc57abfbaced59cef208e6f5130793f7d176b91",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "f9569fafbbbb9e8a7ccf07b7974fa5457b147e3acef989b22f728c4ad63c33fa",
    },
    # Open-board binding and polarization results extend the evaluation.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "79fd7900f8b7dd8a8ba6aac02dc57abfbaced59cef208e6f5130793f7d176b91",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "d66e869b9473071496bc7ec07f5a7ef6392941862d66935fced07a5a222b203a",
    },
    # Common-reference constituent measurements extend the evaluation.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "79fd7900f8b7dd8a8ba6aac02dc57abfbaced59cef208e6f5130793f7d176b91",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "4654816e559c321f11c85942662ede04272bd971ebf07ee4f60a6d0497b5b0a6",
    },
    # Fixed-volume and collar measurements extend the evaluation.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "79fd7900f8b7dd8a8ba6aac02dc57abfbaced59cef208e6f5130793f7d176b91",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "019f7debc2dd8743744d492e312fae9cb36c2b46cad10319b44c8c23b11b081f",
    },
    {
        "model/candidate_laws/interval_echo/SPEC.md": "79fd7900f8b7dd8a8ba6aac02dc57abfbaced59cef208e6f5130793f7d176b91",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "207c1d9e1916b53efbe40549e8590f0ce4c6d13625365783a190bf6f4c756c69",
    },
    # Correlation-energy and concurrent evidence updates replay the complete native message bank.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "79fd7900f8b7dd8a8ba6aac02dc57abfbaced59cef208e6f5130793f7d176b91",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "93a7bbe5c5f3d87a31d2d88fbe2ae8f5dee0d7f87feba76124bd340f9d1d3ba0",
    },
    # The signed area and independent clock observations extend the reviewed evaluation.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "79fd7900f8b7dd8a8ba6aac02dc57abfbaced59cef208e6f5130793f7d176b91",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "acddc0f6282bca3ff37745d040669dc61e277d2a0f16616105fc82d955ba526b",
    },
    # Published candidate descriptions and scorecard, with the native message bank replayed.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "79fd7900f8b7dd8a8ba6aac02dc57abfbaced59cef208e6f5130793f7d176b91",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "7d9f6735d7dd0e59bfcd8bd23354e7cc6398fe805717a6a6cc7044537e05af0b",
    },
    # Guarded-repair scope and same-source Einstein tests update the prose
    # and scorecard. The law and the registered message instrument still replay.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "04e1f8f53881087f159833ed674e38fbbd883f5b5f765a100d369c3ce30d9ee0",
    },
    # Independent wall and dispersion measurements add scorecard evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "ef1b15dd6da778d677998dcb857c8a43c68ff83138bf015f3f53ae5e9a7f71f6",
    },
    # Address interventions, clock bounds, pair recoil and repeated pointer reads extend the evaluation.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "36c4fbcc8c04cf433653a025c0813482064755334a31a1e6a9b4cbead8b9bd8a",
    },
    # Local signal deadlines, material records and field probes extend the scorecard.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "693362b0a40866643a4df5430ad77fec9232959bd07e4b6118b13e3f4ef9c3c7",
    },
    # A total-two sweep, receiver deadlines and a local probe extend the evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "c5de332ed7ae0fca11f99aa9dca1d8979043ad10db0c17d3e02b8827a9294452",
    },
    # Fixed material pointers, circulation-sensitive scattering and orbital transition evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "e74e1da93e64f42f10491686e319f960f469956570218fcc885a2c9426577647",
    },
    # Passive echo calibration, the neutral gap clock and pair-plane controls extend the evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "84ce4a51ef898de57868022f545a0fa62084e338b8661a8d1c303554615e9d34",
    },
    # Pair transport and the source-sensitive clock comparison extend the evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "eeac2a3b1769bb5aa1b614b3ec8974c80806a2872de11f70a0d4e7a8ff40d5a3",
    },
    # Twelve-address records, the released trimer and filtered spatial turns extend the evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "2f7f52c24d9e75726da4f0dfdba4e6aba7ed40337d8ebd5103fe6cbba6ce4a0b",
    },
    # Pair recoil, junction transport and the two-aperture return extend the evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "47e38d1e3640dc535c66279da23a0bd1f3fcd9549dcc3c836913f5e83479e68c",
    },
    # Native daughter modes, regional charge readback and eight material instruments extend the evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "1b0e88d6f0241a1d2ca59c1a9795874e4e9ba1b387ab8f56c02212514004ed67",
    },
    # Same native law: neutral profiles, timed scattering, and stationary geometry.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "313e0df745f5022f8d7d75de1648a33042d3aa955cb733c456b975128a7fe780",
    },
    # Fresh seam interference, released wall receivers, joint clocks, and held geometry.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "9b7ba88b3568b23bddf6476a69041e72bcdf5e01cc86e011c279d4c91c5acef9",
    },
    # Adaptive seam services, whole-band charge and native orbital readbacks extend the evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "d2d42dd38b56bf378560a0cdb379b7ee3001f692176d4f2614c622446e897fff",
    },
    # Source swapping, static parity work and the degenerate clock extend the evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "3cf41343d64c96a35a663f1ac8d391bfd67c77f3fa3863296898cb9f358641e1",
    },
    # Census normalization, inventory effects and the bound arm singlet.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "3df50722895ff6cd9610e264fa725a41a4e9e1268653927335b5b779f3c1c174",
    },
    # The premise rows move to docs/registers/OBSERVABLES.json; only cited register paths change.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "92be169ee8dd0a86a2461204a602e7c7df24151b53d031e6395469da79559e65",
    },
    # The observables ledger review renames premise rows and merges duplicates.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "ee83471e89798c9d7d210e6af90c42127f33cc7ce281ff98eabd52482984f40f",
    },
    # Thermal and screen premises, detector bounds and the strand impulse evidence.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "82045c375d5aef8cb522847e903e8fd452a6cad6fe79adc9447493d2b0f539c7",
    },
    # Variational inputs and the two static interaction measurements.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "72721fcf2747be736d5e2db0d76e19a42d492aa3034c36ca63e113b533fdd569",
    },
    # Platform-independent sampling regenerates receipts that the evaluation cites.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "e660a3752851def2d26862306c78b83e65942756ca662e4443fc6a72101872ec",
    },
    # Charged-vacancy detector, spatial rotation control and free-spectrum bounds.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "138e99ca4157b157389fc89fe1c9009ffdaed25bdde76c3796638ee03aa573ba",
    },
    # Charged formation-energy, magnetic rotations and coherent neutral conversion.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "c0c0a3409d3c254b591c01eb1e73fe0cd00426f3630e0faff2dcf087758fa967",
    },
    # Ground energies, the large parent spectrum and exact register ratios.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "caeac93a6a6015136d705bb9e938ce4a19d784590321cbceba12cf4c56e70b5a",
    },
    # Conjugate recoil, moving trimer, compact parents and resonant delivery.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "38760c952ff6154241f0de6a7e2ec55be640bf24a50c51411480df9586902311",
    },
    # Charged-band rotation, exchange, bound-core recoil and neutral emission.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "6e3ceef7a9be8c35e7fef5f5ffa7ef1879388da0c3db2fcc9fb3a81f8fb7f8dd",
    },
    # Fixed-reference modular and optical response, release ef7c54d1a.
    {
        "model/candidate_laws/interval_echo/SPEC.md": "283ba8edc3a924b12fc71faf54c7daa5a4b90e919892a446c13fde268978808d",
        "model/candidate_laws/interval_echo/candidate.py": "ab125475fd254ce3900a3f71643f9081de37fafe626632788c5e3c1bb5e61106",
        "model/candidate_laws/interval_echo/EVALUATION.json": "dcf774210fc350529bdfd1a128bfc92323044a28c6fc3b83ced3bc04610954f1",
    },
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sources():
    frozen = SPEC["sources"]
    actual = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in frozen}
    reviewed = frozen | REPLAY_UPDATES
    radial = reviewed | REPLAY_70E58E24
    regional = radial | REPLAY_REGIONAL_CONSENSUS
    standard = regional | REPLAY_MODULAR_STANDARD
    gravity = standard | REPLAY_GRAVITY_CAP
    nested = gravity | REPLAY_NESTED_CAP
    joint = nested | REPLAY_JOINT_READ
    directional = joint | REPLAY_DIRECTIONAL_CHARGE
    if actual in (frozen, reviewed, radial, regional, standard, gravity, nested, joint, directional):
        return
    if any(actual == (directional | update) for update in REPLAY_EVALUATIONS):
        return
    mismatches = [name for name, expected in reviewed.items() if actual.get(name) != expected]
    raise ValueError(f"unverified replay source: {', '.join(mismatches)}")


def encode(board, station):
    contact = next(edge for edge in board.external if edge[0] < 12)
    _, a = contact
    b = min(p for p in board.neighbours[a] if p // 12 == a // 12 and board.degree[p] == 5)
    launch = min(set(board.neighbours[a]) & set(board.neighbours[b]))
    require(station["cells"] == board.cells and list(contact) == station["contact"]
            and [a, b] == station["ports"] and launch == station["launch"]
            and launch // 12 == a // 12 != 0, "station map")
    law = echo.native_law
    # Address one seam with common start |1,0> and the same contact area.
    pulse = expm(-1j * np.pi * law.seam_term(echo.KAPPA))[:, law.SRC.index((1, 0))]
    columns = np.zeros((board.n, 2), complex)
    columns[launch, :] = pulse[law.SRC.index((1, 0))]
    columns[[a, b], [0, 1]] = pulse[law.SRC.index((0, 1))]
    require(abs(abs(columns[a, 0]) - 1) < TOL and
            np.linalg.norm(np.delete(pulse, law.SRC.index((0, 1)))) < TOL, "sender exchange")
    return columns


def validate_law(p):
    require(p.shape == (2, 13) and p.dtype.kind == "f" and np.isfinite(p).all()
            and p.min() >= 0 and np.max(abs(p.sum(axis=1) - 1)) <= TOL,
            "incomplete local law")


def readout(h, launch, s):
    state = expm_multiply((-1j * s) * h, launch) if s else launch
    mass = abs(state)**2
    # Audit the raw complement and total before the inherited instrument's guard.
    residual = mass.sum(axis=0) - 1
    outside = mass[12:].sum(axis=0)
    require(np.max(abs(residual)) <= TOL and (1 - mass[:12].sum(axis=0)).min() >= 0,
            "raw normalization or complement")
    p = local_tables(state, np.arange(12)[None, :])[0].T
    validate_law(p)
    require(np.max(abs(outside - p[:, -1])) <= TOL, "outside probability")
    return p, residual.tolist()


def distance(p):
    validate_law(p)
    return float(abs(p[0] - p[1]).sum() / 2)

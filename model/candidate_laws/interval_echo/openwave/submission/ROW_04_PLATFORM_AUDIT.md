# Independent audit of the M10.3 planning and result commits

**Final result head `22458e7f91b394628cd92a5db4eb677e8719d7bc` passes scientific review.** The official M10 mark is 🚧. The row and package are on HOLD for the release conditions in the [completed RER checklist](CHECKLIST.md). The final attachment digests and replays appear in the last section.

**Review on 21 November 2026 UTC.** The scientific reproduction passes. A local platform planning commit, `a0907a62b39af24eb72cd0e2c2101c09520f0e7a`, places the numeric gates and all inputs before the in-platform reproduction. Its immediate parent is the public base `26fbb5fdd07049625d6bb7930af2f6bf545fb6fc`. The task identifies the local release and all results as previously known. This follows the chronology permitted for an existing author result in M10.2. [Official roadmap](https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/m10_roadmap.md#L43) · [M10.2 planning and reproduction record](https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/tasks/m10_2_task_details.md)

I exported the five scientific inputs directly from the platform planning commit into an empty directory. The checkout's sparse index made a repository-wide archive unsuitable, so I used the paths' Git blob objects. I independently reconstructed the author pin in another clean local clone from the public parent and the staged Git bundle. Each of its eleven files matched the source manifest. The author's additional checker compared every matrix entry of both 28,680-state sector Hamiltonians and verified both terminal branches for all four local binary inputs.

| Run on the platform input export | Result | Elapsed |
| --- | --- | ---: |
| `row_04_check.py --board FEDERATION_LEVEL0.npz` | Exit 0. All 24,576 simulated outcomes reproduced, including native counts 4,048, 4,052 and 4,028 and control counts 1,780, 1,182 and 824. | 9.41 s |
| `row_04_window.py --board FEDERATION_LEVEL0.npz` | Exit 0. Minimum mesh probability 0.978883326. Continuous probability bound of at least 94%. | 6.63 s |
| `row_04_check.py --board FEDERATION_LEVEL0.npz --remove-contrast` | Exit 1 at the native persistence gate. | 6.56 s |

The local input freeze, tested law, external board, controls and reported ratios agree. I separately decoded all zlib-compressed single-run streams and compared each of their histograms with its saved count array. The source bundle restores to `4fe10f9ac3e81423155fb3ac7baaf0fdb4c0b11c`. Its eleven additions under `openwave/` leave the native law and geometry identical to the public law commit `7a11a741dcc8de00a2ca826b3800909566467d62`. The archived supplemental encounter at the same cell is additional author evidence and was outside the platform reproduction gates. [Author result](ROW_04.md) · [Author freeze](ROW_04_FREEZE.json)

### Filing checks encountered

1. The planning version placed the registered task in a new table after Conventions. Use the official **In Progress** table in the final result diff so the column's live index finds the task. The result commit can correct the index without changing any earlier gate.
2. OpenWave requires a truthful DCO sign-off on each submitted commit. The reviewed planning commit has no sign-off trailer. The contributor must supply the sign-off for the actual filing branch before it is sent, preserving the registration-before-replay evidence. The auditor cannot grant another person's certification or release authority. [Public contribution instructions](https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/CONTRIBUTING.md#sign-off--developer-certificate-of-origin-dco)
3. The active nomination is a single candidate for the column, with only row 4 requested for scoring. The code and original data on the earlier M10.2 pin remain historical. The proposed source version requires the owner's nomination, coordinated publication, and anonymous source verification before upstream submission. The optional new release using a future upstream commit as a selector has not been drawn.

These are the checks of the planning version. The final package verdict binds any corrections, the exact result commit, the outgoing PR body and the attachment hashes.

## First result packet reviewed

Result head `6ca8851b2f90caaaf6d73352e9fc9e16713a054a` had the planning commit as its immediate parent. I read the two-commit diff, task, canonical candidate entry, roadmap, method, proposed cell, timing receipt and outgoing cover. The result fixed the roadmap index and requested one pin for M10 and full credit only on row 4. It preserved the existing official marks pending maintainer review. The three platform executions agreed with the independent export above. A fresh `git format-patch --full-index --binary --stdout` reproduced the two-commit patch byte-for-byte. The patch passed the Git whitespace checker.

| File in version first reviewed | SHA-256 |
| --- | --- |
| M10_3_PR.md | `b65ac52faaee3babe5afcebc175e96935993842421b829c2a8a50a78b5100ab7` |
| ROW_04_PLATFORM.patch | `1ba3535d0021eeb1191d52dcfe9d10437284196503af03e65d498cd59d96e387` |
| ROW_04_PLATFORM.bundle | `826b6651594ff0a7546a77991b2e556e1ea4915d870228a9e50991d8227bc40f` |
| ROW_04_SOURCE.bundle | `301a24977bfe5d8a7390b994d14c8d5d0ceecb056c7f5c38316971b03238d633` |
| SOURCE_SHA256SUMS | `787f6da6ff90e554826f42391801a7309557f033524531e395b61d94c99e9e0d` |

**Scientific review: PASS at full written criterion.** The exact result and method support the free persistence of this excitation. The gates and rejection have concrete scopes, and the requested scientific decision matches the evidence.

**First distribution review: HOLD.** The method accurately credited the author and identified the unmodified author fixture and two receipt files as CC BY-NC-SA 4.0. The upstream project states that commercial use is permitted and its DCO requires appropriate contribution rights. The named data licence restricts commercial use. Attribution cannot give broader permission to reuse the three files. A practical remedy keeps those files in the pinned author repository under its existing terms and supplies their paths to the in-platform code. M10.2 describes its own use of author data under a pinned checkout. The independent auditor requested a revised clean branch preserving its own pre-replay planning commit. [Source terms](https://github.com/FloatingPragma/observer-patch-holography/blob/7a11a741dcc8de00a2ca826b3800909566467d62/LICENSE) · [OpenWave contribution terms](https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/CONTRIBUTING.md) · [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

Both first-version commits lacked DCO sign-off. The owner decides nomination and publication and the actual contributor signs for the version sent. Source pin `4fe10f9a…` is inspectable through the verified local bundle. Its advertised immutable URL must supply the same files after publication and before filing.

## Earlier revision using separately pinned author data

**Independent scientific and source-distribution review: PASS.** Planning commit `9ebb10fcda45d3769d635de7cd023dc55025a322` follows the public base directly. Result commit `c4fde45b07b1facf649d23c62a92a4479ad4f267` follows that planning commit. The eight changed paths include three Apache-licensed replay scripts, the task and method, roadmap and canonical entries, and a newly generated platform execution record. Neither commit redistributes the RER fixture or either archived author receipt. The platform task reads those inputs offline from the separately pinned author checkout and checks all three digests. Source attribution acknowledges the replay software and points to the author's split licensing. The earlier distribution finding is resolved.

I exported the three Python files directly from the new planning commit into another empty directory and used the independently reconstructed author checkout. The input digests passed. Running the prescribed commands reproduced every address and probability (exit 0, 8.86 s), reproduced the 94% continuous lower bound (exit 0, 4.86 s), and rejected contrast removal at the persistence gate (exit 1, 4.49 s). I inspected the new task, numeric gates, received outputs, method, requested 58-word cell, candidate declaration, In Progress entry, exact outgoing cover and source manifest. The cell has 59 whitespace-delimited tokens including its score symbol and fits the 65-word allowance. The result is a genuine reproduction of previously disclosed measurements under terms committed before the new platform run. The optional future held source is outside this request. An independent full-index binary patch export reproduced the exact outgoing patch bytes, and the whitespace check passed.

| Material at the preceding review | SHA-256 |
| --- | --- |
| [Outgoing PR body](M10_3_PR.md) | `054840b2d30bab285e06be52c8dff10f1726bd87f69af5b3cd3b3ba0c647e807` |
| [Two-commit platform patch](ROW_04_PLATFORM.patch) | `b84eeec4ccebc01cd47b23be2036e9a1e7788b88eb80e7dd3fa398e87a38a96b` |
| [Platform bundle](ROW_04_PLATFORM.bundle) | `a533fc7e9e0d40b98b771d5d3b07be37fd61229696662ab312456b9af35ed0ec` |
| [RER source bundle](ROW_04_SOURCE.bundle) | `301a24977bfe5d8a7390b994d14c8d5d0ceecb056c7f5c38316971b03238d633` |
| [Eleven-file RER manifest](SOURCE_SHA256SUMS) | `787f6da6ff90e554826f42391801a7309557f033524531e395b61d94c99e9e0d` |

**Release disposition.** Full scientific credit is warranted for the named prepared-pair stability criterion. Owner nomination for the precise source, release authorization, truthful DCO sign-off by the actual contributor and public source verification remain pending. Both locally prepared platform commits are unsigned. If signing creates new commit IDs, commit the signed design before rerunning, update the planning references and rebind the exact result, patch and cover. The final cover identifies this audit, the [component checklist](ROW_04_REVIEW.md), method and receipts as its accompanying review attachment. Maintainers decide the official cell.

## Final RER-compliant layout and retest

At the final review the reference public OpenWave base was `26fbb5fdd07049625d6bb7930af2f6bf545fb6fc`. New planning commit `52a7cfe62adbb9009737b502c65d3c12508aa1c3` follows it directly. The final result commit is `22458e7f91b394628cd92a5db4eb677e8719d7bc`. The same staged RER source is now also named by the bundle tag `interval-echo-m10.3-source`. I inspected all eight changed paths, compared the three execution-script blobs to the independently run copies and checked the final full-index binary format-patch against the actual two-commit history. The scripts reside under the registered research scripts directory. The canonical declaration identifies the candidate state, address, free law, memory, preparation, observer, calibration, source ancestry and exclusions. The method now maps its equations to immutable source locations. The integrator's two platform linters passed.

The exact final release script reproduced the counts and historical probabilities in 10.90 s, exit 0. The second script reproduced the continuous bound in 6.54 s, exit 0. The actual removal of the binding interaction rejected at the scientific gate in 6.21 s, exit 1. I also replayed the supplemental 20- and 320-cell comparison using its own source and receipt. The 320-cell full-board probability uncertainty of at most 0.002157368 separates native 0.99098250 from control 0.45000173 at exposure 1. That experiment is a research size check. The platform's 94% window bound applies to the registered 20-cell release.

| Final checked material | SHA-256 |
| --- | --- |
| [Outgoing PR body](M10_3_PR.md) | `be766eeb7897584116cb30561c211e7238f09f2391c11bfdd69f023d21b056f5` |
| [Platform patch](ROW_04_PLATFORM.patch) | `b6241bf4b1b8b942528869aa24371aee4529966b5eab95e131915e94cbd10e96` |
| [Platform bundle](ROW_04_PLATFORM.bundle) | `14895bf76b83fa81f79226a78084641c57546dd7b606512d8d9e549bdf81d5b0` |
| [Tagged RER source bundle](ROW_04_SOURCE.bundle) | `1d3ba519dc7e0e2156ec42784dc98f2d8bcdeba211cb76c970e4d9ec4e6cd4ba` |
| [RER manifest](SOURCE_SHA256SUMS) | `787f6da6ff90e554826f42391801a7309557f033524531e395b61d94c99e9e0d` |

**Scientific G6 and public-source distribution review: PASS for the full written finite-persistence criterion. Formal row and package disposition: HOLD.** The RER guidelines issued at public commit `962ca3d087afa63a7eb8d5c1814e74ffe41f285a` require owner nomination and pushed registration evidence for the release replay. The local planning chain is available. No public pre-run timestamp, owner authorization or DCO sign-off for this version has been recorded. The owner can authorize publication of the signed design first, then the short disclosed-result reproduction, updated version binding, public source verification and PR. [RER checklist and rechecks](CHECKLIST.md)

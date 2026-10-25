# Scientific Methods & Mathematical Reference — MALDI-GlycoID

## 1. Monoisotopic Mass Calculation
- **Released Native Glycan Neutral Mass:**
  $$M_{\text{released}} = \sum_{i} n_i \cdot m_{\text{residue}, i} + m_{\text{H}_2\text{O}} + \Delta m_{\text{label}}$$
  where $m_{\text{H}_2\text{O}} = 18.0105646865\text{ Da}$.
- **Peptide-Attached Glycan Moiety Mass (for Glycopeptides):**
  $$M_{\text{attached}} = \sum_{i} n_i \cdot m_{\text{residue}, i}$$
  omitting the terminal reducing-end $\text{H}_2\text{O}$ lost upon N-glycosidic linkage to Asn.
- **Permethylation:** Calculated by exact site-counting elemental substitution ($\text{OH}/\text{NH} \rightarrow \text{OCH}_3/\text{NCH}_3$, $+14.0156500642\text{ Da}$ per methylated site) for both free reducing-end (`PERMETHYLATED`) and reduced alditol (`PERMETHYLATED_REDUCED`) glycans.
- **Linkage-Specific Sialic Acid Derivatization (`EEC_SIALIC`, `DMA_SIALIC`):** Distinguishes $\alpha 2,3$-linked sialic acids (intramolecular lactonization, $-18.010565\text{ Da}$ or amidation) from $\alpha 2,6$-linked sialic acids (ethyl esterification $+28.031300\text{ Da}$ or dimethylamidation $+27.047285\text{ Da}$).

## 2. Ion `m/z` & Mass Error Formulas
- **Theoretical `m/z`:**
  $$m/z = \frac{M_{\text{neutral}} + \sum c_k \cdot m_{\text{carrier}, k}}{|z|}$$
  where each charged carrier mass $m_{\text{carrier}, k}$ already accounts for the electron mass ($0.000548579909\text{ Da}$).
- **Signed & Absolute Mass Error:**
  $$\Delta m_{\text{Da}} = m/z_{\text{obs}} - m/z_{\text{theo}}, \qquad \Delta m_{\text{ppm}} = \frac{m/z_{\text{obs}} - m/z_{\text{theo}}}{m/z_{\text{theo}}} \times 10^6$$

## 3. Heuristic Ranking Score vs. Isotope Similarity (`idotp`)
- Candidate ranking uses an explicit, transparent heuristic score on $[0, 100]$ combining mass accuracy, relative peak intensity, and optional isotope envelope dot-product (`idotp`).
- This heuristic score is **never** presented as a calibrated posterior probability or false discovery rate.

## 4. Relative Quantification & Missingness
- **Base-Peak Relative Intensity (%):** $100 \times I_{\text{peak}} / I_{\text{max}}$.
- **Relative Glycan Abundance (%):** $100 \times I_{\text{glycan}} / \sum_{g} I_{g}$.
- Missing glycan peaks across replicates are preserved as `None` (never silently zero-filled).

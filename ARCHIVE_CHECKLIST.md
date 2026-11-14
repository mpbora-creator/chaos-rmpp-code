# GitHub–Zenodo release checklist

1. Confirm that all coauthors approve the MIT license and creator order.
2. Replace `OWNER/REPOSITORY` in `README.md`, `CITATION.cff`, and
   `CODE_AVAILABILITY.md` with the final GitHub path.
3. Add ORCID identifiers to `CITATION.cff` and `.zenodo.json` if available.
4. Run `./scripts/audit_release.sh` and `./scripts/reproduce_figures.sh`.
5. Create a new public GitHub repository and push this directory.
6. Sign in to Zenodo, enable the repository under GitHub integration, and
   optionally reserve a DOI before publication.
7. Create an annotated Git tag and GitHub release named `v1.0.0`.
8. Confirm that Zenodo ingests the release and that title, creators, license,
   description, keywords, version, and files are correct.
9. Insert the version-specific Zenodo DOI into the manuscript and
   `CODE_AVAILABILITY.md`.
10. If raw checkpoints or Vlasov snapshots are deposited separately, link the
    software and dataset records using Zenodo's “is supplemented by” relation.

Suggested commands after creating the empty GitHub repository:

```bash
git init
git add .
git commit -m "Archive release v1.0.0"
git branch -M main
git remote add origin https://github.com/OWNER/REPOSITORY.git
git push -u origin main
git tag -a v1.0.0 -m "Version 1.0.0"
git push origin v1.0.0
```

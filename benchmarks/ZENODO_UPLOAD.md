# Zenodo upload of the Kimonu records

Everything is in `@OUT@/`, built from Kimonu.py commit @COMMIT@ by `benchmarks/zenodo_bundle.sh`.

| file | what to do with it |
|---|---|
| `upload/Kimonu_records.zip` (@ZIPMB@ MB) | upload |
| `upload/README.md` | upload (Zenodo shows it as a preview, and it is the same file as `README.md` inside the zip) |
| `zenodo_metadata.json` | the metadata in the format of Zenodo's API, to copy field by field into the form |
| `description.html` | the description, open it in a browser, select all and copy it into the Description field |
| `bundle/` | the unpacked archive, for a look before uploading (do not upload) |
| `work/logs/bundle_*.log` | the build log, which ends with the rebuild of the figures, tables and numbers from the zip |

The record has 2 files and @ZIPMB@ MB, well within Zenodo's limits of 100 files and 50 GB per record.

## Steps

1. Copy the two files of `upload/` to your PC, for example with WinSCP, or in PowerShell with
   `scp pschneider@<cluster login node>:@OUT@/upload/* .`
2. Open https://zenodo.org, log in, and click **New upload** (the **+** menu at the top right).
3. **Files.** Drag `Kimonu_records.zip` and `README.md` into the upload area and wait until both show 100%.
4. **Digital Object Identifier.** At "Do you already have a DOI for this upload?" choose **No**, then click
   **Get a DOI now!**. Zenodo reserves a DOI of the form 10.5281/zenodo.NNNNNNN. Note it down. It resolves only after
   step 15.
5. **Resource type.** Dataset.
6. **Title.** @TITLE@
7. **Publication date.** The day you publish (the form proposes today).
8. **Creators.** Click **Add creator**, choose **Person**, family name `Schneider`, given names `Philipp C.`, ORCID
   `0000-0001-8858-428X`, affiliation `Massachusetts General Hospital` (pick the entry from the list). The form takes
   the institution only, and the department is in `CITATION.cff` and in `zenodo_metadata.json`.
9. **Description.** Paste from `description.html`. It names the article DOI as 10.1126/sciadv.XXXXXXX until the
   article has one, so replace it then or delete those words before publishing.
10. **Licenses.** Creative Commons Attribution 4.0 International (remove any other license the form proposes).
11. **Keywords and subjects.** Add one by one @KEYWORDS@.
12. **Languages** English. **Version** for example 1.0 (`zenodo_metadata.json` has a placeholder).
13. **Related works.** Add
    - relation **Is supplement to**, identifier `https://github.com/VonAlphaBisZulu/Kimonu.py`, scheme URL, resource
      type Software,
    - relation **Is compiled by**, identifier `https://github.com/VonAlphaBisZulu/Kimonu.py/tree/@COMMIT_FULL@`,
      scheme URL, resource type Software (the commit must be pushed to GitHub for this link to work),
    - once the article has its DOI, relation **Is supplement to**, the article DOI, scheme DOI, resource type
      Publication / Journal article. The metadata can be edited after publishing without a new version, so this row
      can wait.
14. **References.** Paste the three lines of `references` in `zenodo_metadata.json`.
15. **Visibility** Public. Click **Save draft**, then **Preview** and check the page, then **Publish**. After publishing
    the files cannot be changed (a correction needs **New version**), but the metadata can.

Science Advances asks that the data be available when the paper is published. Publishing the record at submission is
the simplest way to let the reviewers open it.

## After the DOI exists

- Manuscript (`Kimonu.md`, data availability), replace "==(DOI)==" by `(https://doi.org/10.5281/zenodo.NNNNNNN)`, and
  fill in the release tag where the Materials and Methods ask for the version number.
- Cover letter, the same DOI wherever the archive is mentioned.
- Optional, before step 15. The README inside the zip carries placeholders for the archive DOI and the article DOI,
  and CITATION.cff also for the version and date of this record. The README names Kimonu by commit and points to the
  releases of the repository for the latest stable version. To fill them, rebuild after step 4 (about 30 minutes, on the
  cluster), delete the two uploaded files in the draft and upload the new ones from the same place.

      cd ~/pipelines/2026-08_Kimonu/Kimonu.py
      KIMONU_REV=@COMMIT@ ZENODO_DOI=10.5281/zenodo.NNNNNNN ARCHIVE_VERSION=1.0 RELEASE_DATE=YYYY-MM-DD \
        bash benchmarks/zenodo_bundle.sh @OUT@

  `ARTICLE_DOI=10.1126/sciadv...` can be added the same way once it exists. The same command
  without the variables rebuilds the whole package after any change of the figures, Data S1 and S2 or the records.

-- ---------------------------------------------------------------------------
-- Build a run-level sample sheet from the public SRA metadata table in
-- BigQuery (nih-sra-datastore.sra.metadata).
--
-- Column reference:
--   https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-metadata-table/
--
-- Sample attributes (tissue, sex, treatment, genotype, ...) are submitter-
-- defined and live in a nested key/value array, so they need UNNEST rather
-- than a fixed column. Attribute keys are NOT standardised across studies:
-- run the audit query at the bottom first to see what a given study actually
-- used, then edit the CASE block below to match.
--
-- Parameters are declared for the bq CLI / python client:
--   bq query --use_legacy_sql=false --parameter=study:STRING:SRP123456 \
--     < resources/sra_metadata.sql
-- ---------------------------------------------------------------------------

DECLARE study STRING DEFAULT @study;

WITH runs AS (
  SELECT
    m.acc                AS run,
    m.sample_name        AS sample_name,
    m.biosample          AS biosample,
    m.sra_study          AS study,
    m.bioproject         AS bioproject,
    m.organism           AS organism,
    m.assay_type         AS assay_type,
    m.librarylayout      AS library_layout,
    m.libraryselection   AS library_selection,
    m.platform           AS platform,
    m.instrument         AS instrument,
    m.avgspotlen         AS avg_read_length,
    m.mbases             AS mbases,
    -- Flatten the attribute array into a lookup we can pull named keys from.
    (SELECT ARRAY_AGG(STRUCT(LOWER(a.k) AS k, a.v AS v)) FROM UNNEST(m.attributes) AS a) AS attrs
  FROM `nih-sra-datastore.sra.metadata` AS m
  WHERE m.sra_study = study
    AND m.consent = 'public'
    AND m.assay_type = 'RNA-Seq'
),

annotated AS (
  SELECT
    run,
    sample_name,
    biosample,
    study,
    bioproject,
    organism,
    library_layout,
    platform,
    instrument,
    avg_read_length,
    mbases,
    -- COALESCE across the aliases submitters commonly use for the same field.
    COALESCE(
      (SELECT v FROM UNNEST(attrs) WHERE k = 'treatment_sam'  LIMIT 1),
      (SELECT v FROM UNNEST(attrs) WHERE k = 'condition_sam'  LIMIT 1),
      (SELECT v FROM UNNEST(attrs) WHERE k = 'genotype_sam'   LIMIT 1),
      'UNSPECIFIED'
    ) AS condition,
    COALESCE(
      (SELECT v FROM UNNEST(attrs) WHERE k = 'tissue_sam'      LIMIT 1),
      (SELECT v FROM UNNEST(attrs) WHERE k = 'source_name_sam' LIMIT 1)
    ) AS tissue,
    COALESCE(
      (SELECT v FROM UNNEST(attrs) WHERE k = 'sex_calc' LIMIT 1),
      (SELECT v FROM UNNEST(attrs) WHERE k = 'sex_sam'  LIMIT 1)
    ) AS sex
  FROM runs
)

SELECT
  run,
  -- Snakemake wildcards choke on spaces and slashes; normalise sample labels.
  REGEXP_REPLACE(COALESCE(sample_name, run), r'[^A-Za-z0-9_.-]', '_') AS sample,
  REGEXP_REPLACE(LOWER(condition), r'[^a-z0-9_.-]', '_')             AS condition,
  library_layout,
  tissue,
  sex,
  organism,
  platform,
  instrument,
  avg_read_length,
  mbases,
  biosample,
  study,
  bioproject
FROM annotated
ORDER BY condition, sample;


-- ---------------------------------------------------------------------------
-- AUDIT QUERY — run this first on a new study to see which attribute keys
-- exist and how many runs carry each one. Attribute naming is inconsistent
-- across submitters, and guessing is how sample sheets end up silently wrong.
--
--   SELECT LOWER(a.k) AS attribute_key,
--          COUNT(DISTINCT m.acc) AS n_runs,
--          ARRAY_AGG(DISTINCT a.v LIMIT 10) AS example_values
--   FROM `nih-sra-datastore.sra.metadata` AS m, UNNEST(m.attributes) AS a
--   WHERE m.sra_study = 'SRP123456'
--   GROUP BY attribute_key
--   ORDER BY n_runs DESC;
-- ---------------------------------------------------------------------------

--
-- MaSim database creation script, as generated by pgAdmin
--

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

CREATE SCHEMA sim;

ALTER SCHEMA sim OWNER TO sim;

SET default_tablespace = '';

SET default_with_oids = false;

CREATE TABLE sim.monthlydata (
    id integer NOT NULL,
    replicateid integer NOT NULL,
    dayselapsed integer NOT NULL,
    modeltime bigint NOT NULL,
    seasonalfactor integer NOT NULL,
    treatmentfailures double precision NOT NULL,
    beta double precision NOT NULL,
    entrytime timestamp with time zone
);


ALTER TABLE sim.monthlydata OWNER TO sim;

CREATE TABLE sim.monthlysitedata (
    monthlydataid integer NOT NULL,
    locationid integer NOT NULL,
    population integer NOT NULL,
    clinicalepisodes integer NOT NULL,
    treatments integer NOT NULL,
    treatmentfailures integer NOT NULL,
    eir double precision NOT NULL,
    pfprunder5 double precision NOT NULL,
    pfpr2to10 double precision NOT NULL,
    pfprall double precision NOT NULL
);


ALTER TABLE sim.monthlysitedata OWNER TO sim;

CREATE VIEW public.v_runningstats AS
 SELECT running.id,
    running.dayselapsed,
    (running.dayselapsed / 365) AS years,
    running.popuation,
    running.seasonalfactor,
    (running.endtime - running.starttime) AS "time"
   FROM ( SELECT md.id,
            md.seasonalfactor,
            md.dayselapsed,
            md.entrytime AS endtime,
            ( SELECT monthlydata.entrytime
                   FROM sim.monthlydata
                  WHERE (monthlydata.id = (md.id - 1))) AS starttime,
            ( SELECT sum(monthlysitedata.population) AS population
                   FROM sim.monthlysitedata
                  WHERE (monthlysitedata.monthlydataid = md.id)) AS popuation
           FROM sim.monthlydata md) running
  ORDER BY running.id;


ALTER TABLE public.v_runningstats OWNER TO postgres;

CREATE TABLE sim.configuration (
    id integer NOT NULL,
    yaml character varying NOT NULL,
    md5 character varying NOT NULL,
    name character varying,
    notes character varying,
    filename character varying,
    ncols integer default -1 NOT NULL,
    nrows integer default -1 NOT NULL,
    xllcorner double precision default 0 NOT NULL,
    yllcorner double precision default 0 NOT NULL,
    cellsize integer default -1 NOT NULL
);


ALTER TABLE sim.configuration OWNER TO sim;

CREATE SEQUENCE sim.configuration_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE sim.configuration_id_seq OWNER TO sim;

ALTER SEQUENCE sim.configuration_id_seq OWNED BY sim.configuration.id;

CREATE TABLE sim.genotype (
    id integer NOT NULL,
    name character varying NOT NULL
);


ALTER TABLE sim.genotype OWNER TO sim;

CREATE TABLE sim.location (
    id integer NOT NULL,
    configurationid integer NOT NULL,
    index integer NOT NULL,
    x integer NOT NULL,
    y integer NOT NULL,
    beta double precision NOT NULL,
    district integer
);


ALTER TABLE sim.location OWNER TO sim;

CREATE SEQUENCE sim.location_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE sim.location_id_seq OWNER TO sim;

ALTER SEQUENCE sim.location_id_seq OWNED BY sim.location.id;

CREATE SEQUENCE sim.monthlydata_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE sim.monthlydata_id_seq OWNER TO sim;

ALTER SEQUENCE sim.monthlydata_id_seq OWNED BY sim.monthlydata.id;

CREATE TABLE sim.monthlygenomedata (
    monthlydataid integer NOT NULL,
    locationid integer NOT NULL,
    genomeid integer NOT NULL,
	occurrences integer NOT NULL,
    clinicaloccurrences integer NOT NULL,
    occurrences0to5 integer NOT NULL,
    occurrences2to10 integer NOT NULL,
    weightedfrequency double precision NOT NULL
);


ALTER TABLE sim.monthlygenomedata OWNER TO sim;

CREATE TABLE sim.replicate (
    id integer NOT NULL,
    configurationid integer NOT NULL,
    seed bigint NOT NULL,
    starttime timestamp with time zone NOT NULL,
    endtime timestamp with time zone,
    movement character(1)
);


ALTER TABLE sim.replicate OWNER TO sim;

CREATE SEQUENCE sim.replicate_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE sim.replicate_id_seq OWNER TO sim;

ALTER SEQUENCE sim.replicate_id_seq OWNED BY sim.replicate.id;

CREATE TABLE sim.study (
    id integer NOT NULL,
    name character varying NOT NULL
);


ALTER TABLE sim.study OWNER TO sim;

CREATE SEQUENCE sim.study_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE sim.study_id_seq OWNER TO sim;

ALTER SEQUENCE sim.study_id_seq OWNED BY sim.study.id;

CREATE TABLE sim.studyconfigurations (
    studyid integer NOT NULL,
    configurationid integer NOT NULL
);

ALTER TABLE sim.studyconfigurations OWNER TO sim;

ALTER TABLE ONLY sim.configuration ALTER COLUMN id SET DEFAULT nextval('sim.configuration_id_seq'::regclass);

ALTER TABLE ONLY sim.location ALTER COLUMN id SET DEFAULT nextval('sim.location_id_seq'::regclass);

ALTER TABLE ONLY sim.monthlydata ALTER COLUMN id SET DEFAULT nextval('sim.monthlydata_id_seq'::regclass);

ALTER TABLE ONLY sim.replicate ALTER COLUMN id SET DEFAULT nextval('sim.replicate_id_seq'::regclass);

ALTER TABLE ONLY sim.study ALTER COLUMN id SET DEFAULT nextval('sim.study_id_seq'::regclass);

ALTER TABLE ONLY sim.configuration
    ADD CONSTRAINT configuration_pkey PRIMARY KEY (id);

ALTER TABLE ONLY sim.genotype
    ADD CONSTRAINT genotype_pkey PRIMARY KEY (id);

ALTER TABLE ONLY sim.location
    ADD CONSTRAINT location_index_unique UNIQUE (id, configurationid, index);

ALTER TABLE ONLY sim.location
    ADD CONSTRAINT location_pkey PRIMARY KEY (id);

ALTER TABLE ONLY sim.monthlydata
    ADD CONSTRAINT monthlydata_pkey PRIMARY KEY (id);

ALTER TABLE ONLY sim.monthlygenomedata
    ADD CONSTRAINT monthlygenomedata_pkey PRIMARY KEY (monthlydataid, genomeid, locationid);

ALTER TABLE ONLY sim.monthlysitedata
    ADD CONSTRAINT monthlysitedata_pkey PRIMARY KEY (monthlydataid, locationid);

ALTER TABLE ONLY sim.replicate
    ADD CONSTRAINT replicate_pkey PRIMARY KEY (id);

ALTER TABLE ONLY sim.study
    ADD CONSTRAINT study_pkey PRIMARY KEY (id);

ALTER TABLE ONLY sim.studyconfigurations
    ADD CONSTRAINT studyconfigurations_pkey PRIMARY KEY (studyid, configurationid);

CREATE INDEX fki_g ON sim.monthlygenomedata USING btree (genomeid);

CREATE INDEX fki_location_configurationid_fk ON sim.location USING btree (configurationid);

CREATE INDEX fki_monthlydata_replicateid_fk ON sim.monthlydata USING btree (replicateid);

CREATE INDEX fki_monthlygenomedata_monthlydataid_fk ON sim.monthlygenomedata USING btree (monthlydataid);

CREATE INDEX fki_replicate_studyid_fk ON sim.replicate USING btree (configurationid);

CREATE INDEX fki_studyconfigurations_configuratinid_fk ON sim.studyconfigurations USING btree (configurationid);

CREATE INDEX fki_studyconfigurations_studyid_fk ON sim.studyconfigurations USING btree (studyid);

ALTER TABLE ONLY sim.monthlysitedata
    ADD CONSTRAINT "MonhtlyDataId_FK" FOREIGN KEY (monthlydataid) REFERENCES sim.monthlydata(id);

ALTER TABLE ONLY sim.location
    ADD CONSTRAINT location_configurationid_fk FOREIGN KEY (configurationid) REFERENCES sim.configuration(id);

ALTER TABLE ONLY sim.monthlydata
    ADD CONSTRAINT monthlydata_replicateid_fk FOREIGN KEY (replicateid) REFERENCES sim.replicate(id);

ALTER TABLE ONLY sim.monthlygenomedata
    ADD CONSTRAINT monthlygenomedata_genotypeid_fk FOREIGN KEY (genomeid) REFERENCES sim.genotype(id);

ALTER TABLE ONLY sim.monthlygenomedata
    ADD CONSTRAINT monthlygenomedata_monthlydataid_fk FOREIGN KEY (monthlydataid) REFERENCES sim.monthlydata(id);

ALTER TABLE ONLY sim.replicate
    ADD CONSTRAINT replicate_configurationid_fk FOREIGN KEY (configurationid) REFERENCES sim.configuration(id);

ALTER TABLE ONLY sim.studyconfigurations
    ADD CONSTRAINT studyconfigurations_configuratinid_fk FOREIGN KEY (configurationid) REFERENCES sim.configuration(id);

ALTER TABLE ONLY sim.studyconfigurations
    ADD CONSTRAINT studyconfigurations_studyid_fk FOREIGN KEY (studyid) REFERENCES sim.study(id);

-- The fine movement table isn't strictly necessary for the database
CREATE SEQUENCE sim.movement_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

CREATE TABLE sim.movement
(
    id bigint NOT NULL DEFAULT nextval('sim.movement_id_seq'::regclass),
    replicateid integer NOT NULL,
    timestep integer NOT NULL,
    individualid integer NOT NULL,
    source integer NOT NULL,
    destination integer NOT NULL,
    CONSTRAINT movement_pkey PRIMARY KEY (id),
    CONSTRAINT movement_replicateid_fk FOREIGN KEY (replicateid)
        REFERENCES sim.replicate (id) MATCH SIMPLE
        ON UPDATE NO ACTION
        ON DELETE NO ACTION
        NOT VALID
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE sim.movement OWNER to sim;
ALTER SEQUENCE sim.movement_id_seq OWNER to sim;

CREATE INDEX fki_movement_replicateid_fk
    ON sim.movement USING btree
    (replicateid)
    TABLESPACE pg_default;

-- The coarse movement table isn't strictly necessary for the database
CREATE SEQUENCE sim.districtmovement_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

CREATE TABLE sim.districtmovement
(
    id integer NOT NULL DEFAULT nextval('sim.districtmovement_id_seq'::regclass),
    replicateid integer NOT NULL,
    timestep integer NOT NULL,
    count integer NOT NULL,
    source integer NOT NULL,
    destination integer NOT NULL,
    CONSTRAINT districtmovement_pkey PRIMARY KEY (id),
    CONSTRAINT districtmovement_replicateid_fk FOREIGN KEY (replicateid)
        REFERENCES sim.replicate (id) MATCH SIMPLE
        ON UPDATE NO ACTION
        ON DELETE NO ACTION
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE sim.districtmovement OWNER to sim;
ALTER SEQUENCE sim.districtmovement_id_seq OWNER to sim;

CREATE INDEX fki_districtmovement_replicateid_fk
    ON sim.districtmovement USING btree
    (replicateid)
    TABLESPACE pg_default;
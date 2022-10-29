-- Creating DB
CREATE DATABASE test -- Set the name of your new database here
    WITH
    OWNER = superuser -- Add in your own superuser here
;

-- Use DB
\c test -- This must match the name you used in line 2!

-- Add module for uuid
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Creating tables
CREATE TABLE constant (
	const_id UUID NOT NULL DEFAULT uuid_generate_v1() PRIMARY KEY,
	value NUMERIC,
	precision INT,
	time_added timestamp DEFAULT current_timestamp
);

CREATE TABLE named_constant (
    const_id UUID NOT NULL PRIMARY KEY REFERENCES constant (const_id),
	name VARCHAR NOT NULL UNIQUE,
	description VARCHAR,
	artificial INT NOT NULL DEFAULT 0
);

CREATE TABLE pcf_canonical_constant (
    const_id UUID NOT NULL PRIMARY KEY REFERENCES constant (const_id),
	P INT[] NOT NULL,
	Q INT[] NOT NULL,
	last_matrix TEXT, -- solely because of the absurdly huge numbers that can go here... not even NUMERIC is enough...
	depth INT,
	convergence INT,
	
	UNIQUE(P, Q)
);

CREATE TABLE scan_history (
    const_id UUID NOT NULL PRIMARY KEY REFERENCES constant (const_id),
	algorithm VARCHAR NOT NULL,
	time_scanned timestamp DEFAULT current_timestamp,
	details VARCHAR
);

CREATE TABLE relation (
	relation_id UUID NOT NULL DEFAULT uuid_generate_v1() PRIMARY KEY,
	relation_type VARCHAR NOT NULL,
	details INT[] NOT NULL,
	time_added timestamp DEFAULT current_timestamp
);

CREATE TABLE constant_in_relation (
	const_id UUID NOT NULL REFERENCES constant (const_id) ON UPDATE CASCADE ON DELETE CASCADE,
	relation_id UUID NOT NULL REFERENCES relation (relation_id) ON UPDATE CASCADE,
	CONSTRAINT const_relation_pkey PRIMARY KEY (const_id, relation_id)
);

CREATE TABLE relation_audit (
	relation_id UUID NOT NULL REFERENCES relation (relation_id),
	operation CHAR(1) NOT NULL,
	stamp timestamp NOT NULL DEFAULT current_timestamp,
	userid TEXT NOT NULL
);

-- based on https://www.postgresql.org/docs/current/plpgsql-trigger.html
CREATE OR REPLACE FUNCTION process_relation_audit() RETURNS TRIGGER AS $relation_audit$
	BEGIN
		IF (TP_OG = 'INSERT') THEN
			INSERT INTO relation_audit SELECT NEW.relation_id, 'I', now(), user;
		ELSIF (TP_OG = 'UPDATE') THEN
			INSERT INTO relation_audit SELECT NEW.relation_id, 'U', now(), user;
		ELSIF (TP_OG = 'DELETE') THEN
			INSERT INTO relation_audit SELECT OLD.relation_id, 'D', now(), user;
		END IF;
		RETURN NULL;
	END;
$relation_audit$ LANGUAGE plpgsql;

CREATE TRIGGER relation_audit
AFTER INSERT OR UPDATE OR DELETE ON relation
	FOR EACH ROW EXECUTE FUNCTION process_relation_audit();

DROP ROLE IF EXISTS scout;
CREATE ROLE scout WITH
	NOLOGIN
	NOSUPERUSER
	INHERIT
	NOCREATEDB
	NOCREATEROLE
	NOREPLICATION;

REVOKE ALL ON constant FROM scout;
REVOKE ALL ON constant_in_relation FROM scout;
REVOKE ALL ON named_constant FROM scout;
REVOKE ALL ON pcf_canonical_constant FROM scout;
REVOKE ALL ON relation FROM scout;
REVOKE ALL ON relation_audit FROM scout;
REVOKE ALL ON scan_history FROM scout;

GRANT SELECT ON constant TO scout;
GRANT SELECT ON constant_in_relation TO scout;
GRANT SELECT ON named_constant TO scout;
GRANT SELECT ON pcf_canonical_constant TO scout;
GRANT SELECT ON relation TO scout;
-- no, you're not allowed to see the audits table >:|
GRANT SELECT ON scan_history TO scout;

GRANT REFERENCES ON named_constant TO scout;
GRANT REFERENCES ON pcf_canonical_constant TO scout;

GRANT INSERT ON constant_in_relation TO scout;
GRANT INSERT ON relation TO scout;

-- Then when someone new wants to contribute, run code similar to this:
-- CREATE ROLE [username] LOGIN;
-- ALTER USER [username] WITH PASSWORD '[password]'
-- GRANT scout to [username];
-- also just in case every password works, see https://stackoverflow.com/a/21054627

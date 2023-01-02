-- Creating DB
CREATE DATABASE test -- Set the name of your new database here
    WITH
    OWNER = superuser -- Add in your own superuser here
;

-- Use DB, name here must match the name you used in line 2!
\c test

-- Add module for uuid
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Creating tables
CREATE TABLE constant (
	const_id UUID NOT NULL DEFAULT uuid_generate_v1() PRIMARY KEY,
	value NUMERIC,
	precision INT,
	time_added timestamp DEFAULT current_timestamp
	priority INT NOT NULL DEFAULT 1
	tweeted INT NOT NULL DEFAULT 0
);

CREATE TABLE named_constant (
    const_id UUID NOT NULL PRIMARY KEY REFERENCES constant (const_id),
	name VARCHAR NOT NULL UNIQUE,
	description VARCHAR,
);

CREATE TABLE pcf_canonical_constant (
    const_id UUID NOT NULL PRIMARY KEY REFERENCES constant (const_id),
    original_a INT[],
    original_b INT[],
	"P" INT[] NOT NULL, -- must be in quotes otherwise it becomes lowercase...
	"Q" INT[] NOT NULL,
	last_matrix TEXT, -- solely because of the absurdly huge numbers that can go here... not even NUMERIC is enough...
	depth INT,
	convergence INT,
	
	UNIQUE("P", "Q")
);

CREATE TABLE derived_constant (
    const_id UUID NOT NULL PRIMARY KEY REFERENCES constant (const_id),
	family VARCHAR NOT NULL,
	args JSONB NOT NULL
);

CREATE TABLE pcf_family (
    family_id UUID NOT NULL DEFAULT uuid_generate_v1() PRIMARY KEY,
	a VARCHAR NOT NULL,
	b VARCHAR NOT NULL,
	
	UNIQUE(a, b)
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
	precision INT,
	time_added timestamp DEFAULT current_timestamp
	priority INT NOT NULL DEFAULT 1
	tweeted INT NOT NULL DEFAULT 0
);

CREATE TABLE constant_in_relation (
	const_id UUID NOT NULL REFERENCES constant (const_id) ON UPDATE CASCADE ON DELETE CASCADE,
	relation_id UUID NOT NULL REFERENCES relation (relation_id) ON UPDATE CASCADE ON DELETE CASCADE,
	CONSTRAINT const_relation_pkey PRIMARY KEY (const_id, relation_id)
);

REVOKE ALL ON SCHEMA public FROM PUBLIC; -- public is overprivileged yo

DROP ROLE IF EXISTS spectator;
CREATE ROLE spectator WITH
	NOLOGIN
	NOSUPERUSER
	INHERIT
	NOCREATEDB
	NOCREATEROLE
	NOREPLICATION;

GRANT SELECT ON constant TO spectator;
GRANT SELECT ON constant_in_relation TO spectator;
GRANT SELECT, REFERENCES ON named_constant TO spectator;
GRANT SELECT, REFERENCES ON pcf_canonical_constant TO spectator;
GRANT SELECT, REFERENCES ON derived_constant TO spectator;
GRANT SELECT ON pcf_family TO spectator;
GRANT SELECT ON relation TO spectator;
GRANT SELECT ON scan_history TO spectator;


DROP ROLE IF EXISTS spectator_public;
CREATE ROLE spectator_public WITH
	LOGIN
	NOSUPERUSER
	INHERIT
	NOCREATEDB
	NOCREATEROLE
	NOREPLICATION;

GRANT spectator to spectator_public;
ALTER USER spectator_public WITH PASSWORD 'helloworld123'; -- exploit amazon RDS forbidding changing passwords, using rds.restrict_password_commands parameter


DROP ROLE IF EXISTS scout;
CREATE ROLE scout WITH
	NOLOGIN
	NOSUPERUSER
	INHERIT
	NOCREATEDB
	NOCREATEROLE
	NOREPLICATION;

GRANT spectator TO scout;
GRANT INSERT ON constant_in_relation TO scout;
GRANT INSERT ON relation TO scout;


DROP ROLE IF EXISTS pioneer;
CREATE ROLE pioneer WITH
	NOLOGIN
	NOSUPERUSER
	INHERIT
	NOCREATEDB
	NOCREATEROLE
	NOREPLICATION;

GRANT scout TO pioneer;
GRANT INSERT ON constant TO pioneer;
GRANT INSERT ON pcf_canonical_constant TO pioneer;
GRANT INSERT ON derived_constant TO pioneer;
GRANT INSERT ON pcf_family TO pioneer;


DROP ROLE IF EXISTS janitor;
CREATE ROLE janitor WITH
	NOLOGIN
	NOSUPERUSER
	INHERIT
	NOCREATEDB
	NOCREATEROLE
	NOREPLICATION;

GRANT scout TO janitor;
GRANT UPDATE ON constant TO janitor;
GRANT UPDATE ON pcf_canonical_constant TO janitor;
GRANT UPDATE ON derived_constant TO janitor;
GRANT UPDATE, DELETE ON relation TO janitor;
GRANT UPDATE, DELETE ON constant_in_relation TO janitor;


DROP ROLE IF EXISTS twitterbot;
CREATE ROLE twitterbot WITH
	NOLOGIN
	NOSUPERUSER
	INHERIT
	NOCREATEDB
	NOCREATEROLE
	NOREPLICATION;

GRANT scout TO twitterbot;
GRANT UPDATE (tweeted) ON constant TO twitterbot;
GRANT UPDATE (tweeted) ON relation TO twitterbot;

-- Then when someone new wants to contribute, run code similar to this:
-- CREATE ROLE [username] LOGIN;
-- ALTER USER [username] WITH PASSWORD '[password]'
-- GRANT [scout/pioneer] to [username];
-- also just in case every password works, see https://stackoverflow.com/a/21054627

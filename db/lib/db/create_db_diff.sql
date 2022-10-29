-- Creating DB
CREATE DATABASE ramanujan
    WITH
    OWNER = postgres
;

-- Use DB
\c ramanujan

CREATE TABLE cf_precision (
	cf_id UUID NOT NULL PRIMARY KEY REFERENCES cf (cf_id),
	insertion_date timestamp DEFAULT current_timestamp,
	depth INT NOT NULL,
	precision INT NOT NULL,
	value NUMERIC NOT NULL,
	previous_calc NUMERIC[] NOT NULL,
	general_data JSONB,
	interesting NUMERIC DEFAULT 0,
	update_time timestamp DEFAULT current_timestamp
);

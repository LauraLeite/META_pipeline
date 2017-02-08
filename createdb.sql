create database META_database;

use META_database;

create table gene (CONTIG varchar(255) NOT NULL, START_MATCH int NOT NULL, END_MATCH int NOT NULL, READ_ID varchar(255) NOT NULL  PRIMARY KEY, GENE_ID varchar(255) NOT NULL, GENE_LENGTH int NOT NULL, FOREIGN KEY(GENE_ID) REFERENCES ko(GENE_ID));

load data local infile 'gene.txt' into table amostra_gene;

create table tax (READ_ID varchar(255) NOT NULL PRIMARY KEY, GIDESC varchar(255));

load data local infile 'blastn.txt' into table blastn;

create table uni (UNI_ID varchar(255), KO varchar(255) NOT NULL PRIMARY KEY);

load data local infile 'Uni2KO' into table uni;

create table kodesc (CAT1 varchar(255), CAT2 varchar(255), KO varchar(155) NOT NULL PRIMARY KEY, DESCRIP varchar(555));

load data local infile 'KEGG_hier.txt' into table kodesc;

create table ko (GENE_ID varchar(255) NOT NULL PRIMARY KEY, UNI varchar(50));

load data local infile 'ko.txt' into table ko;

alter table gene add column RATIO float;

alter table ko add column KO varchar;

update ko set uni.UNI_ID = ko.UNI from ko q inner join uni a on q.ko = a.uni where q.KO is null

update gene set RATIO = ABS(END_MATCH-START_MATCH)/GENE_LENGTH;

CREATE VIEW sample_view AS select k.CONTIG,k.START_MATCH,k.END_MATCH,k.READ_ID,k.GENE_ID,k.GENE_LENGTH,k.RATIO,g.KO,kt.GIDESC,gp.CAT1,gp.CAT2,gp.DESCRIP from ko g inner join gene k on k.GENE_ID=g.GENE_ID  inner join blastn kt on kt.READ_ID=k.READ_ID inner join kodesc gp on gp.KO =g. KO ;

select CAT2,GENE_ID,KO,GIDESC,sum(RATIO),count(*) from sample_view group by CAT2,KO,GIDESC;

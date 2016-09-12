-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

-- -----------------------------------------------------
-- Schema mydb
-- -----------------------------------------------------
-- -----------------------------------------------------
-- Schema mapper_dev
-- -----------------------------------------------------

-- -----------------------------------------------------
-- Table `centers`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `centers` ;

CREATE TABLE IF NOT EXISTS `centers` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `name` VARCHAR(45) NOT NULL,
  `created_at` DATETIME NOT NULL,
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `studies`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `studies` ;

CREATE TABLE IF NOT EXISTS `studies` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `name` VARCHAR(45) NOT NULL,
  `created_at` DATETIME NOT NULL,
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `hosts`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `hosts` ;

CREATE TABLE IF NOT EXISTS `hosts` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `name` VARCHAR(45) NOT NULL,
  `created_at` DATETIME NOT NULL,
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `pis`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `pis` ;

CREATE TABLE IF NOT EXISTS `pis` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `name` VARCHAR(45) NOT NULL,
  `created_at` DATETIME NOT NULL,
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `projects`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `projects` ;

CREATE TABLE IF NOT EXISTS `projects` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `name` VARCHAR(45) NULL DEFAULT NULL,
  `created_at` DATETIME NOT NULL,
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `samples`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `samples` ;

CREATE TABLE IF NOT EXISTS `samples` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `sample_id` VARCHAR(45) NOT NULL,
  `center_id` INT(11) NOT NULL,
  `study_id` INT(11) NOT NULL,
  `pi_id` INT(11) NOT NULL,
  `host_id` INT(11) NOT NULL,
  `project_id` INT(11) NOT NULL,
  `filename` VARCHAR(45) NOT NULL,
  `run_dir` VARCHAR(45) NOT NULL,
  `fullpath` TEXT NOT NULL,
  `year` INT NULL,
  `reads` BIGINT NULL,
  `created_at` DATETIME NOT NULL,
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `states`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `states` ;

CREATE TABLE IF NOT EXISTS `states` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `name` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE INDEX `name_UNIQUE` (`name` ASC),
  INDEX `index3` (`name` ASC))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `fastqs`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `fastqs` ;

CREATE TABLE IF NOT EXISTS `fastqs` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `sample_id` INT(11) NOT NULL,
  `build` VARCHAR(45) NOT NULL DEFAULT '38',
  `path` TEXT NOT NULL,
  `read_group` TEXT NOT NULL,
  `created_at` DATETIME NOT NULL,
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `results`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `results` ;

CREATE TABLE IF NOT EXISTS `results` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `sample_id` INT(11) NOT NULL,
  `state_id` INT(11) NOT NULL,
  `build` VARCHAR(45) NOT NULL DEFAULT '38',
  `exported_at` DATETIME NULL DEFAULT NULL,
  `created_at` DATETIME NOT NULL,
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `steps`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `steps` ;

CREATE TABLE IF NOT EXISTS `steps` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `name` VARCHAR(45) NOT NULL DEFAULT 'all',
  PRIMARY KEY (`id`),
  UNIQUE INDEX `name_UNIQUE` (`name` ASC),
  INDEX `index3` (`name` ASC))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `jobs`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `jobs` ;

CREATE TABLE IF NOT EXISTS `jobs` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `result_id` INT(11) NOT NULL,
  `job_id` INT(11) NOT NULL,
  `step_id` INT(11) NOT NULL,
  `cluster` VARCHAR(45) NOT NULL,
  `procs` INT(11) NOT NULL,
  `memory` INT(11) NOT NULL,
  `walltime` VARCHAR(45) NOT NULL,
  `exit_code` INT(11) NULL DEFAULT NULL,
  `elapsed` INT(11) NULL DEFAULT '0',
  `node` VARCHAR(45) NULL DEFAULT NULL,
  `delay` INT(11) NULL DEFAULT '0',
  `submitted_at` DATETIME NULL DEFAULT NULL,
  `started_at` DATETIME NULL DEFAULT NULL,
  `ended_at` DATETIME NULL DEFAULT NULL,
  `created_at` DATETIME NOT NULL,
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  INDEX `index4` (`cluster` ASC),
  INDEX `index5` (`job_id` ASC))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `logs`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `logs` ;

CREATE TABLE IF NOT EXISTS `logs` (
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `job_id` INT(11) NOT NULL,
  `timestamp` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `level` VARCHAR(45) NOT NULL,
  `message` TEXT NOT NULL,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = latin1;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;

-- -----------------------------------------------------
-- Data for table `states`
-- -----------------------------------------------------
START TRANSACTION;
INSERT INTO `states` (`id`, `name`) VALUES (DEFAULT, 'cancelled');
INSERT INTO `states` (`id`, `name`) VALUES (DEFAULT, 'completed');
INSERT INTO `states` (`id`, `name`) VALUES (DEFAULT, 'failed');
INSERT INTO `states` (`id`, `name`) VALUES (DEFAULT, 'requested');
INSERT INTO `states` (`id`, `name`) VALUES (DEFAULT, 'started');
INSERT INTO `states` (`id`, `name`) VALUES (DEFAULT, 'submitted');
INSERT INTO `states` (`id`, `name`) VALUES (DEFAULT, 'converted');

COMMIT;


-- -----------------------------------------------------
-- Data for table `steps`
-- -----------------------------------------------------
START TRANSACTION;
INSERT INTO `steps` (`id`, `name`) VALUES (DEFAULT, 'all');
INSERT INTO `steps` (`id`, `name`) VALUES (DEFAULT, 'bam2fastq');
INSERT INTO `steps` (`id`, `name`) VALUES (DEFAULT, 'align');
INSERT INTO `steps` (`id`, `name`) VALUES (DEFAULT, 'cloud-align');

COMMIT;


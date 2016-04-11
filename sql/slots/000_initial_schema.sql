-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

-- -----------------------------------------------------
-- Schema mydb
-- -----------------------------------------------------
-- -----------------------------------------------------
-- Schema slots_dev
-- -----------------------------------------------------

-- -----------------------------------------------------
-- Table `projects`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `projects` ;

CREATE TABLE IF NOT EXISTS `projects` (
  `id` INT(11) NOT NULL AUTO_INCREMENT COMMENT '',
  `name` VARCHAR(45) NOT NULL COMMENT '',
  `created_at` DATETIME NOT NULL COMMENT '',
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP COMMENT '',
  PRIMARY KEY (`id`)  COMMENT '',
  UNIQUE INDEX `index2` (`name` ASC)  COMMENT '')
ENGINE = InnoDB
AUTO_INCREMENT = 2
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `types`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `types` ;

CREATE TABLE IF NOT EXISTS `types` (
  `id` INT(11) NOT NULL AUTO_INCREMENT COMMENT '',
  `name` VARCHAR(45) NOT NULL COMMENT '',
  `created_at` DATETIME NOT NULL COMMENT '',
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP COMMENT '',
  PRIMARY KEY (`id`)  COMMENT '',
  UNIQUE INDEX `index2` (`name` ASC)  COMMENT '')
ENGINE = InnoDB
AUTO_INCREMENT = 2
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `pools`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `pools` ;

CREATE TABLE IF NOT EXISTS `pools` (
  `id` INT(11) NOT NULL AUTO_INCREMENT COMMENT '',
  `project_id` INT(11) NOT NULL COMMENT '',
  `type_id` INT(11) NOT NULL COMMENT '',
  `name` VARCHAR(45) NOT NULL COMMENT '',
  `hostname` VARCHAR(45) NOT NULL COMMENT '',
  `size_used` BIGINT(20) NOT NULL DEFAULT '0' COMMENT '',
  `size_total` BIGINT(20) NOT NULL DEFAULT '0' COMMENT '',
  `path` VARCHAR(45) NOT NULL COMMENT '',
  `active` TINYINT(4) NOT NULL DEFAULT '1' COMMENT '',
  `created_at` DATETIME NOT NULL COMMENT '',
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP COMMENT '',
  PRIMARY KEY (`id`)  COMMENT '',
  UNIQUE INDEX `index6` (`name` ASC, `hostname` ASC)  COMMENT '',
  INDEX `pools_idx_project_id` (`project_id` ASC)  COMMENT '',
  INDEX `pools_idx_type_id` (`type_id` ASC)  COMMENT '',
  CONSTRAINT `pools_fk_project_id`
    FOREIGN KEY (`project_id`)
    REFERENCES `projects` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `pools_fk_type_id`
    FOREIGN KEY (`type_id`)
    REFERENCES `types` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB
AUTO_INCREMENT = 3
DEFAULT CHARACTER SET = latin1;


-- -----------------------------------------------------
-- Table `slots`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `slots` ;

CREATE TABLE IF NOT EXISTS `slots` (
  `id` INT(11) NOT NULL AUTO_INCREMENT COMMENT '',
  `pool_id` INT(11) NOT NULL COMMENT '',
  `name` VARCHAR(255) NOT NULL COMMENT '',
  `sha1` CHAR(40) NULL DEFAULT NULL COMMENT '',
  `size` VARCHAR(45) NOT NULL COMMENT '',
  `created_at` DATETIME NOT NULL COMMENT '',
  `modified_at` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP COMMENT '',
  PRIMARY KEY (`id`)  COMMENT '',
  UNIQUE INDEX `index3` (`name` ASC, `pool_id` ASC)  COMMENT '',
  INDEX `slots_idx_pool_id` (`pool_id` ASC)  COMMENT '',
  CONSTRAINT `slots_fk_pool_id`
    FOREIGN KEY (`pool_id`)
    REFERENCES `pools` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB
AUTO_INCREMENT = 3
DEFAULT CHARACTER SET = latin1;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;

-- -----------------------------------------------------
-- Data for table `types`
-- -----------------------------------------------------
START TRANSACTION;
INSERT INTO `types` (`id`, `name`, `created_at`, `modified_at`) VALUES (DEFAULT, 'nfs', DEFAULT, DEFAULT);

COMMIT;

